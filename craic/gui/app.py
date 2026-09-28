"""CRAIC main window — wires the canvas, panels, engines, and ambiguity tools."""

from __future__ import annotations

import os
import time
from typing import Optional

import numpy as np
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDialog, QDockWidget,
    QFileDialog, QInputDialog, QLabel, QLineEdit,
    QMainWindow, QMenu, QMessageBox, QProgressDialog, QPushButton, QSlider,
    QTabWidget, QToolBar, QToolButton,
)

from .. import AUTHOR_URL, CITATION, WEBSITE, accel, editing, figures, io, session as session_mod
from ..ambiguity import disagreement as disagree_mod
from ..ambiguity import ensemble as ensemble_mod
from ..ambiguity import posterior as post_mod
from ..ambiguity import reliability as rel_mod
from ..ambiguity import report as report_mod
from ..ambiguity import trimming as trim_mod
from ..domain import (Alignment, Alphabet, CodingSpec, IdMismatch, Level,
                       code_name, genetic_codes)
from ..engines import CodonAware, available_engines, builtin_variants
from . import colors
from .canvas import AlignmentCanvas
from . import tracks as track_mod
from .dialogs import (EngineParamsDialog, FigureExportDialog, JobProgress, SimulateDialog,
                      _fmt_params, duration)
from .document import Document
from .panels import ColumnInspector, PosteriorPanel, SandboxPanel, ScoreLegend
from .workers import run_async, run_cancellable

_LEVEL_LABEL = {Level.NT: "Nucleotide", Level.CODON: "Codon", Level.AA: "Amino acid"}
# The column overlays, described once in craic.gui.tracks — see that module for
# why. These module-level names are the compatibility surface the rest of the
# window (and the tests) use.
_TRACK_DEFS = track_mod.build(lambda aln: _conservation(aln), trim_mod)
_TRACKS = [t.label for t in _TRACK_DEFS]
_TRUTH_TRACK = track_mod.index_of(_TRACK_DEFS, "truth")
_CHEAP_TRACKS = frozenset(i for i, t in enumerate(_TRACK_DEFS) if t.cheap and i != 0)
# The residue-colouring modes, in the order the combo lists them.
_COLOUR_RESIDUE, _COLOUR_CONFIDENCE, _COLOUR_TRUTH = 0, 1, 2
# Said wherever a long analysis starts without the compiled core.
_NO_CORE = ("CRAIC's compiled core is not installed, so this runs tens of times more "
            "slowly than it should. Reinstalling CRAIC usually fixes that.")


def _close_job_dialog(job: dict) -> None:
    """Close a job's progress dialog, once, without it reporting a cancel:
    closing a QProgressDialog emits ``canceled``."""
    dlg = job.pop("dialog", None)
    if dlg is None:
        return
    dlg.canceled.disconnect()
    dlg.close()
    dlg.deleteLater()


def _conservation(aln) -> np.ndarray:
    """Per-column conservation: the fraction of *all* sequences that share the
    commonest residue. Gaps count in the denominator (a gap is a state of its own —
    a 'fifth base' / 21st residue), so a column that is 80% T and 20% gap scores
    0.8, not 1.0. All-gap columns carry no signal and are left as NaN."""
    from collections import Counter

    n = aln.n_seqs or 1
    out = np.full(aln.length, np.nan)
    for c in range(aln.length):
        residues = [r[c].upper() for r in aln.rows if r[c] not in "-.?"]
        if residues:
            out[c] = Counter(residues).most_common(1)[0][1] / n
    return out


def _friendly_open_error(path, exc):
    """Turn a parser exception (which may embed an entire offending sequence) into a
    short, readable popup message: returns (title, text)."""
    import re

    name = os.path.basename(path)
    raw = " ".join(str(exc).split())            # collapse newlines / runs of whitespace
    if "illegal character" in raw.lower():
        m = re.search(r"[Ii]llegal character[: ]+'?(\S)'?", raw)
        ch = m.group(1) if m else "?"
        t = re.search(r"[Tt]axon (\S+?):", raw)
        who = f'Sequence "{t.group(1)}"' if t else "A sequence"
        return (f"Can't open {name}",
                f'{who} contains the character "{ch}", which is not a nucleotide, '
                "RNA, or protein residue.\n\n"
                "If this is a binary gene-content / presence-absence matrix (0/1) or "
                "another discrete-character dataset, CRAIC can't display it \u2014 it reads "
                "molecular sequence alignments.")
    return f"Can't open {name}", raw[:300] + (" \u2026" if len(raw) > 300 else "")


class CraicWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CRAIC — Conserved-Region Alignment by Iterative Convergence")
        self.resize(1180, 720)

        # What the user is working on — alignment, reference, accuracy, groups,
        # annotations, paths, unsaved state, generation — lives in its own object
        # (see craic.gui.document). This window owns widgets.
        self.doc = Document(self)
        self.reliability: Optional[rel_mod.Reliability] = None
        self.agreement: Optional[np.ndarray] = None
        self._job: Optional[dict] = None            # the long analysis running, if any
        self._reliability_waiters: list = []        # what is waiting for the reliability report
        self._perturbation_skipped = False          # the user settled for consistency alone
        self._consistency_secs: Optional[float] = None
        # Coalesce a burst of edits into one recovery write.
        self._autosave_timer = QTimer(self)
        self._autosave_timer.setSingleShot(True)
        self._autosave_timer.setInterval(4000)
        self._autosave_timer.timeout.connect(self._write_autosave)
        self._align_note = ""
        self._align_progress = None
        self._align_token = None
        self._align_cancelled = False
        self._last_probe = None
        self._hotspot_idx = -1
        self._warned_frame_breaks = False
        self._undo_stack = []
        self._redo_stack = []
        self._engine_params: dict = {}     # engine.key -> {param_key: value}
        self._inspect_col: Optional[int] = None    # column the inspector is describing
        self._truth_view: Optional[dict] = None    # view state parked while the answer is up
        self._engines = available_engines()

        self.canvas = AlignmentCanvas()
        self.canvas.selectionChanged.connect(self._on_selection)
        self.canvas.residueClicked.connect(self._on_residue_clicked)
        self.canvas.editRequested.connect(self._on_edit_key)
        self.canvas.cursorMoved.connect(self._on_cursor_moved)
        self.canvas.rowsSelected.connect(self._on_rows_selected)
        self.canvas.selectionCleared.connect(self._on_deselect)
        self.setCentralWidget(self.canvas)

        self.sandbox = SandboxPanel()
        self.sandbox.params_for = self._params_for     # the ⚙ settings reach the sandbox
        self.sandbox.applied.connect(lambda a: self._set_alignment(a))
        self.posterior = PosteriorPanel()
        self.inspector = ColumnInspector()
        self.inspector.jumpRequested.connect(self._jump_to_column)
        tabs = QTabWidget()
        self.tabs = tabs
        tabs.addTab(self.sandbox, "Realignment sandbox")
        tabs.addTab(self.posterior, "Posterior explorer")
        tabs.addTab(self.inspector, "Column inspector")
        self.dock = QDockWidget("Ambiguity tools", self)
        self.dock.setWidget(tabs)
        self.dock.setMinimumWidth(360)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dock)
        self.dock.hide()   # hidden by default; View menu > "Ambiguity tools panel" shows it

        self._build_toolbars()
        self._build_statusbar()
        self._refresh_controls()

    # ------------------------------------------------------------------ #
    def _build_toolbars(self):
        tb = QToolBar("File")
        self.addToolBar(tb)
        op = QAction("Open…", self); op.setShortcut("Ctrl+O"); op.triggered.connect(self._open)
        sv = QAction("Save alignment as…", self)
        sv.setShortcut("Ctrl+Shift+S"); sv.triggered.connect(self._save)
        # Sessions keep the whole working state; exporting an alignment keeps the
        # residues. Ctrl+S is the session, because that is the one that loses
        # nothing and so is the safe thing for a reflex to reach for.
        sess_save = QAction("Save session", self)
        sess_save.setShortcut("Ctrl+S")
        sess_save.triggered.connect(lambda: self._save_session())
        sess_save_as = QAction("Save session as…", self)
        sess_save_as.triggered.connect(lambda: self._save_session(ask=True))
        sess_open = QAction("Open session…", self)
        sess_open.triggered.connect(self._open_session)
        self._session_actions = [sess_save, sess_save_as]
        ei = QAction("Export image…", self); ei.triggered.connect(self._export_image)
        em = QAction("Export masked…", self); em.triggered.connect(self._export_masked)
        hi = QAction("History…", self); hi.triggered.connect(self._show_history)
        sa = QAction("Select all", self); sa.setShortcut("Ctrl+A")
        sa.triggered.connect(lambda: self.canvas.select_all())
        srt = QAction("Sort by similarity", self); srt.triggered.connect(self._sort_by_similarity)
        deg = QAction("Remove gap columns", self); deg.triggered.connect(self._remove_gap_columns)
        # Copy → one item per writable format; Ctrl+C is the FASTA shortcut.
        copy_menu = QMenu("Copy", self)
        self._copy_actions = []
        for key, label, _ext in io.write_formats():
            a = QAction(label, self)
            a.triggered.connect(lambda _checked=False, k=key: self._copy_as(k))
            if key == "fasta":
                a.setShortcut("Ctrl+C")
            copy_menu.addAction(a)
            self._copy_actions.append(a)
        # Publication figures — one dialog to choose the figure and its inputs.
        fig_act = QAction("Export figure…", self)
        fig_act.triggered.connect(self._export_fig)
        self._fig_actions = [fig_act]
        hot = QAction("Go to next ambiguous region", self)
        hot.setShortcut("Ctrl+G")
        hot.triggered.connect(self._next_hotspot)
        # Manual editing — undo/redo + structural edits, all on the provenance log.
        undo = QAction("Undo", self); undo.setShortcut("Ctrl+Z"); undo.triggered.connect(self._undo)
        redo = QAction("Redo", self); redo.setShortcut("Ctrl+Shift+Z"); redo.triggered.connect(self._redo)
        nl = QAction("Nudge residue left", self); nl.setShortcut("Ctrl+,")
        nl.triggered.connect(lambda: self._edit_nudge(-1))
        nr = QAction("Nudge residue right", self); nr.setShortcut("Ctrl+.")
        nr.triggered.connect(lambda: self._edit_nudge(1))
        ig = QAction("Insert gap column", self); ig.triggered.connect(self._edit_insert_gap_col)
        dc = QAction("Delete empty column", self); dc.triggered.connect(self._edit_delete_col)
        pin = QAction("Pin selected columns", self); pin.setShortcut("Ctrl+P")
        pin.setToolTip("Mark the selected columns as trusted, so a realign keeps them fixed.")
        pin.triggered.connect(self._pin_selected)
        unpin = QAction("Unpin selected columns", self); unpin.setShortcut("Ctrl+Shift+P")
        unpin.triggered.connect(self._unpin_selected)
        unpin_all = QAction("Clear all pins", self); unpin_all.triggered.connect(self._clear_pins)
        ar = QAction("Realign around pinned columns…", self); ar.triggered.connect(self._anchored_realign)
        rm = QAction("Mask selected residues", self); rm.setShortcut("Ctrl+K")
        rm.setToolTip("Mark residues you judge misaligned. The column is kept; on "
                      "export they are written as missing data (N or X).")
        rm.triggered.connect(lambda: self._mask_selected_residues(True))
        ru = QAction("Unmask selected residues", self); ru.setShortcut("Ctrl+Shift+K")
        ru.triggered.connect(lambda: self._mask_selected_residues(False))
        rt = QAction("Mask residues below the threshold", self)
        rt.setToolTip("Mask every residue whose reliability is below the Mask "
                      "slider's value, keeping all columns.")
        rt.triggered.connect(self._mask_residues_below)
        rc = QAction("Clear residue masks", self); rc.triggered.connect(self._clear_residue_masks)
        edit_align = QMenu("Edit alignment", self)
        for a in (undo, redo):
            edit_align.addAction(a)
        edit_align.addSeparator()
        for a in (nl, nr, ig, dc):
            edit_align.addAction(a)
        edit_align.addSeparator()
        for a in (pin, unpin, unpin_all, ar):
            edit_align.addAction(a)
        edit_align.addSeparator()
        for a in (rm, ru, rt, rc):
            edit_align.addAction(a)
        self._edit_actions = [undo, redo, nl, nr, ig, dc, pin, unpin, unpin_all, ar,
                              rm, ru, rt, rc]
        cmp = QAction("Compare with alignment…", self)
        cmp.triggered.connect(self._compare_alignment)
        # Teaching: a dataset whose answer is known, and the answer itself.
        gen = QAction("Generate dataset with known answer…", self)
        gen.triggered.connect(self._generate_dataset)
        ref_open = QAction("Load reference alignment…", self)
        ref_open.triggered.connect(self._load_reference)
        ref_clear = QAction("Forget reference alignment", self)
        ref_clear.triggered.connect(self._clear_reference)
        ref_show = QAction("Show the true alignment", self)
        ref_show.setCheckable(True)
        ref_show.setShortcut("Ctrl+T")
        ref_show.toggled.connect(self._toggle_true_alignment)
        self.show_truth_act = ref_show
        self._truth_actions = [ref_open, ref_clear, ref_show]
        # Keep shortcuts live at the window level even when the menus are closed.
        for a in (op, sv, ei, em, hi, sa, srt, deg, hot, cmp, gen,
                  sess_open, sess_save, sess_save_as,
                  *self._copy_actions, *self._fig_actions, *self._edit_actions,
                  *self._truth_actions):
            self.addAction(a)
        file_menu = QMenu(self)
        file_menu.addAction(op); file_menu.addAction(sess_open); file_menu.addSeparator()
        file_menu.addAction(sess_save); file_menu.addAction(sess_save_as)
        file_menu.addAction(sv); file_menu.addSeparator()
        file_menu.addAction(ei); file_menu.addAction(em); file_menu.addSeparator()
        file_menu.addAction(fig_act); file_menu.addSeparator()
        file_menu.addAction(cmp); file_menu.addSeparator()
        file_menu.addAction(hi)
        # Named sequence groups + column annotations — rebuilt when the menu opens.
        self.groups_menu = QMenu("Sequence groups", self)
        self.groups_menu.aboutToShow.connect(self._rebuild_groups_menu)
        self.annot_menu = QMenu("Annotations", self)
        self.annot_menu.aboutToShow.connect(self._rebuild_annot_menu)
        edit_menu = QMenu(self)
        edit_menu.addMenu(copy_menu); edit_menu.addAction(sa); edit_menu.addSeparator()
        edit_menu.addAction(srt); edit_menu.addAction(deg); edit_menu.addSeparator()
        edit_menu.addAction(hot); edit_menu.addMenu(edit_align)
        edit_menu.addMenu(self.groups_menu); edit_menu.addMenu(self.annot_menu)
        # Teach: everything that needs a known answer to be meaningful.
        teach_menu = QMenu(self)
        teach_menu.addAction(gen)
        teach_menu.addSeparator()
        teach_menu.addAction(ref_open)
        teach_menu.addAction(ref_clear)
        teach_menu.addSeparator()
        teach_menu.addAction(ref_show)
        file_btn = QToolButton()
        file_btn.setText("File")
        file_btn.setPopupMode(QToolButton.InstantPopup)
        file_btn.setMenu(file_menu)
        edit_btn = QToolButton()
        edit_btn.setText("Edit")
        edit_btn.setPopupMode(QToolButton.InstantPopup)
        edit_btn.setMenu(edit_menu)
        # View menu — the Ambiguity tools dock is hidden by default; this toggle
        # shows/hides it (and reflects its current state).
        view_menu = QMenu(self)
        toggle = self.dock.toggleViewAction()
        toggle.setText("Ambiguity tools panel")
        view_menu.addAction(toggle)
        view_menu.addSeparator()
        # Standard zoom shortcuts: Qt maps "Ctrl" to ⌘ on macOS and Ctrl on
        # Windows/Linux, and the menu shows the right symbol for each platform.
        bigger = QAction("Bigger text", self)
        bigger.setShortcuts([QKeySequence.ZoomIn, QKeySequence("Ctrl+=")])
        bigger.triggered.connect(lambda: self.canvas.zoom(1))
        smaller = QAction("Smaller text", self)
        smaller.setShortcut(QKeySequence.ZoomOut)
        smaller.triggered.connect(lambda: self.canvas.zoom(-1))
        reset_txt = QAction("Reset text size", self)
        reset_txt.setShortcut(QKeySequence("Ctrl+0"))
        reset_txt.triggered.connect(self.canvas.reset_zoom)
        self._zoom_actions = [bigger, smaller, reset_txt]
        for a in self._zoom_actions:
            view_menu.addAction(a)
        for a in self._zoom_actions:
            self.addAction(a)        # keep zoom shortcuts live with the menu closed
        view_btn = QToolButton()
        view_btn.setText("View")
        view_btn.setPopupMode(QToolButton.InstantPopup)
        view_btn.setMenu(view_menu)
        teach_btn = QToolButton()
        teach_btn.setText("Teach")
        teach_btn.setPopupMode(QToolButton.InstantPopup)
        teach_btn.setMenu(teach_menu)
        # Help: the website, and who made CRAIC and how to cite it.
        help_menu = QMenu(self)
        web_act = QAction("CRAIC website", self)
        web_act.triggered.connect(lambda: _open_url(WEBSITE))
        about_act = QAction("About CRAIC and how to cite it…", self)
        about_act.triggered.connect(self._show_about)
        help_menu.addAction(web_act); help_menu.addAction(about_act)
        help_btn = QToolButton()
        help_btn.setText("Help")
        help_btn.setPopupMode(QToolButton.InstantPopup)
        help_btn.setMenu(help_menu)
        # Tool buttons default to a smaller font than the neighbouring labels and
        # buttons; match the app font so the whole ribbon reads at one size.
        for b in (file_btn, edit_btn, view_btn, teach_btn, help_btn):
            b.setFont(QApplication.font())
        tb.addWidget(file_btn); tb.addWidget(edit_btn); tb.addWidget(view_btn)
        tb.addWidget(teach_btn); tb.addWidget(help_btn); tb.addSeparator()
        # Actions that only make sense once an alignment is loaded.
        self._doc_actions = [sv, *self._session_actions, ei, em, hi, sa, srt, deg, hot, cmp,
                             *self._copy_actions, *self._fig_actions, *self._edit_actions]

        tb.addWidget(QLabel(" Engine: "))
        self.engine_combo = QComboBox()
        for e in self._engines:
            self.engine_combo.addItem(e.label, e)
        self.engine_combo.setCurrentIndex(self._default_engine_index())
        tb.addWidget(self.engine_combo)
        gear = QPushButton("⚙")
        gear.setToolTip("Edit this aligner's parameters")
        gear.setFixedWidth(30)
        gear.clicked.connect(self._edit_engine_params)
        tb.addWidget(gear)
        self.codon_chk = QCheckBox("align as protein")
        self.codon_chk.setToolTip(
            "Translate the coding sequences, align the amino acids, then "
            "thread the original nucleotides back through that alignment. The result "
            "keeps the nucleotides in memory and is shown as protein; switch the View "
            "level to see codon or nucleotide at any time.")
        # Ticking it changes which alphabet the engine is handed, so which
        # engines can be used changes with it.
        self.codon_chk.toggled.connect(lambda *_: self._sync_engine_combo())
        tb.addWidget(self.codon_chk)
        al = QPushButton("Align"); al.clicked.connect(self._align); tb.addWidget(al)
        self.align_btn = al

        self.addToolBarBreak()
        tb2 = QToolBar("View / Analyze")
        self.addToolBar(tb2)
        tb2.addWidget(QLabel(" View: "))
        self.level_combo = QComboBox()
        # This combo is filled per alignment (Nucleotide / Codon / Amino acid),
        # so at construction it is empty and Qt sizes it to nothing; without
        # these it stays too narrow to show its own entries. AdjustToContents
        # re-sizes on every repopulation, and the minimum contents length keeps
        # it stable rather than jumping between alignments.
        self.level_combo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.level_combo.setMinimumContentsLength(max(len(s) for s in _LEVEL_LABEL.values()))
        self.level_combo.currentIndexChanged.connect(self._on_level)
        tb2.addWidget(self.level_combo)
        self.coding_chk = QCheckBox("treat as coding")
        self.coding_chk.stateChanged.connect(self._on_coding)
        tb2.addWidget(self.coding_chk)
        # The genetic code is not a detail: mitochondrial and Mollicute genes read
        # TGA as tryptophan, and under the standard code every one of those
        # translates to a stop, which the emission model treats as an unknown
        # residue carrying no homology signal.
        self.code_combo = QComboBox()
        for tid, name in genetic_codes():
            self.code_combo.addItem(f"{tid}  {name}", tid)
        self.code_combo.setToolTip(
            "NCBI genetic code, used by the codon and amino-acid views and by "
            "“align as protein”. Mitochondrial and Mollicute (Mycoplasma, "
            "Spiroplasma) sequences do not use the standard code, and reading "
            "them with it turns real tryptophans into stop codons.")
        self.code_combo.currentIndexChanged.connect(self._on_code)
        tb2.addWidget(self.code_combo)
        tb2.addWidget(QLabel(" Colour: "))
        self.color_combo = QComboBox()
        self.color_combo.addItems(["Residue", "Confidence", "Correct placement"])
        self.color_combo.setItemData(
            _COLOUR_TRUTH,
            "Colour every residue by whether it is where the known-true alignment "
            "puts it. Needs a reference — simulate a dataset, or open one.",
            Qt.ToolTipRole)
        self.color_combo.currentIndexChanged.connect(self._on_colour)
        tb2.addWidget(self.color_combo)
        # Colour key for the red-amber-green reliability ramp. Shown whenever
        # something on screen is painted with it — the confidence colouring or a
        # score track — and hidden otherwise so it isn't decoration.
        self.legend = ScoreLegend()
        self.legend.setVisible(False)
        tb2.addWidget(self.legend)
        tb2.addSeparator()

        tb2.addWidget(QLabel(" Track: "))
        self.track_combo = QComboBox(); self.track_combo.addItems(_TRACKS)
        self.track_combo.currentIndexChanged.connect(self._on_track)
        self.track_combo.setToolTip("Overlay a per-column statistic on the alignment. "
                                    "Reliability and Aligner agreement are computed on demand "
                                    "the first time you pick them, then cached.")
        tb2.addWidget(self.track_combo)
        tb2.addSeparator()

        tb2.addWidget(QLabel(" Mask < "))
        self.thr = QSlider(Qt.Horizontal); self.thr.setRange(0, 100); self.thr.setValue(50)
        self.thr.setFixedWidth(120); self.thr.valueChanged.connect(self._on_threshold)
        tb2.addWidget(self.thr)
        self.thr_lbl = QLabel("0.50"); tb2.addWidget(self.thr_lbl)
        self.outlier_chk = QCheckBox("outliers")
        self.outlier_chk.setToolTip("Flag whole sequences that align poorly with the "
                                    "rest (EvalMSA / OD-seq-style); flagged rows get a "
                                    "red marker.")
        self.outlier_chk.stateChanged.connect(self._on_outliers)
        tb2.addWidget(self.outlier_chk)

        self.addToolBarBreak()
        tb3 = QToolBar("Display")
        self.addToolBar(tb3)
        self.scheme_lbl = QLabel(" Scheme: ")
        tb3.addWidget(self.scheme_lbl)
        self.scheme_combo = QComboBox()
        self.scheme_combo.addItems(colors.SCHEME_LABELS)
        self.scheme_combo.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.scheme_combo.currentIndexChanged.connect(self._on_scheme)
        tb3.addWidget(self.scheme_combo)
        self.consensus_chk = QCheckBox("consensus")
        self.consensus_chk.stateChanged.connect(self._on_consensus)
        tb3.addWidget(self.consensus_chk)
        tb3.addSeparator()
        tb3.addWidget(QLabel(" Find: "))
        self.find_edit = QLineEdit()
        self.find_edit.setFixedWidth(160)
        self.find_edit.setPlaceholderText("name or motif")
        self.find_edit.returnPressed.connect(self._find)
        tb3.addWidget(self.find_edit)

    def _build_statusbar(self):
        self.backend_lbl = QLabel()
        self.sel_lbl = QLabel("no selection")
        self.statusBar().addWidget(self.backend_lbl, 1)
        self.statusBar().addPermanentWidget(self.sel_lbl)
        eng = ", ".join(e.label.split(" ")[0] for e in self._engines)
        self.backend_lbl.setText(
            f"core: {accel.backend()}   |   engines: {eng or 'built-in only'}"
        )

    # ------------------------------------------------------------------ #
    def open_path(self, path: str):
        # A session document opens as a session, whatever the Open dialog was
        # called: the user asked for that file, not for a parse of its JSON.
        if session_mod.is_session_file(path):
            self.open_session_path(path)
            return
        try:
            aln = io.load_alignment(path)
        except Exception as exc:
            title, text = _friendly_open_error(path, exc)
            QMessageBox.critical(self, title, text)
            return
        if aln.meta.get("unaligned"):
            QMessageBox.information(
                self, "Unaligned input",
                "Sequences are not aligned. Pick an engine and press Align.")
        aln.meta["history"] = [f"opened {os.path.basename(path)}"]
        # Named sequence groups ride alongside in the sidecar; hold them at the app
        # level (keyed by sequence id) so they survive edits, which rebuild the meta.
        self._groups = dict(aln.meta.get("groups", {}))
        self._annotations = list(aln.meta.get("annotations", []))
        self._path = path
        self._session_path = ""
        self.reference = self.truth = None
        self._reference_label = ""
        self._set_alignment(aln, dirty=False)

    def _open(self):
        if not self._confirm_discard("Opening another file"):
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open sequences / alignment", "",
            "Sequences and sessions (*.fa *.fasta *.fna *.faa *.phy *.aln *.sto *.nex "
            f"*{session_mod.SUFFIX});;All files (*)")
        if path:
            self.open_path(path)

    def _save(self):
        if not self.aln:
            return
        fmts = io.write_formats()  # (key, label, ext)
        filters = ";;".join(f"{label} (*{ext})" for _, label, ext in fmts)
        path, selected = QFileDialog.getSaveFileName(
            self, "Save alignment", "alignment.fasta", filters)
        if not path:
            return
        key, ext = fmts[0][0], fmts[0][2]
        for k, label, e in fmts:
            if selected == f"{label} (*{e})":
                key, ext = k, e
                break
        if not os.path.splitext(path)[1]:
            path += ext
        try:
            self.aln.meta["groups"] = self._groups        # written to the sidecar by io
            self.aln.meta["annotations"] = self._annotations
            io.save_alignment(path, self.aln, key)
            self._path = path
            # The residues, groups and annotations are now on disk (the latter
            # two in the sidecar), so the work is not at risk and the document is
            # clean. What a standard format cannot carry is the truth-mode
            # reference and the view state, so say so rather than let the user
            # discover it tomorrow.
            self._set_dirty(False)
            self._clear_autosave()
            if self.reference is not None:
                self.sel_lbl.setText(
                    f"saved {os.path.basename(path)} — the reference alignment is "
                    "not part of a standard format; use Save session to keep it")
            hist = self.aln.meta.get("history", [])
            if hist:
                side = os.path.splitext(path)[0] + ".craic-history.txt"
                with open(side, "w") as fh:
                    fh.write("CRAIC alignment provenance\n")
                    for i, h in enumerate(hist):
                        fh.write(f"{i + 1}. {h}\n")
        except Exception as exc:
            QMessageBox.critical(self, "Save failed", str(exc))

    def _copy_as(self, fmt: str):
        if not self.aln:
            return
        records = self.canvas.selected_records()
        if not records:
            return
        level = self.canvas.level
        alph = Alphabet.PROTEIN if level == Level.AA else self.aln.alphabet
        sub = Alignment.from_records(records, alphabet=alph)
        try:
            text = io.dump_alignment(sub, fmt)
        except Exception as exc:
            QMessageBox.critical(self, "Copy failed", str(exc)[:400])
            return
        QApplication.clipboard().setText(text)
        label = next(l for k, l, _ in io.write_formats() if k == fmt)
        self.sel_lbl.setText(f"copied selection to clipboard ({label})")

    def _export_image(self):
        if not self.aln:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export alignment image",
                                              "alignment.png", "PNG image (*.png)")
        if not path:
            return
        if not os.path.splitext(path)[1]:
            path += ".png"
        self.canvas.render_full().save(path)

    def _show_history(self):
        if not self.aln:
            return
        hist = self.aln.meta.get("history", [])
        body = "\n".join(f"{i + 1}.  {h}" for i, h in enumerate(hist)) or "No edits recorded yet."
        box = QMessageBox(self)
        box.setWindowTitle("Alignment provenance")
        box.setText("What produced the current alignment, in order\n(saved alongside on Save):")
        box.setInformativeText(body)
        box.setIcon(QMessageBox.Information)
        box.exec()

    def _on_scheme(self, *_):
        self.canvas.set_scheme(self.scheme_combo.currentText().lower())

    def _on_consensus(self, *_):
        self.canvas.set_consensus(self.consensus_chk.isChecked())

    def _sort_by_similarity(self):
        if not self.aln or self.aln.n_seqs < 3:
            return
        from ..progressive import similarity_order
        order = similarity_order([r.replace("-", "") for r in self.aln.rows])
        ids = [self.aln.ids[i] for i in order]
        rows = [self.aln.rows[i] for i in order]
        self._commit(Alignment(ids, rows, self.aln.alphabet, coding=self.aln.coding),
                     "sorted sequences by similarity")

    def _find(self):
        if not self.aln:
            return
        if not self.canvas.search(self.find_edit.text()):
            self.sel_lbl.setText(f"no match for '{self.find_edit.text()}'")

    def _remove_gap_columns(self):
        if not self.aln:
            return
        aln = self.aln
        if aln.coding is not None and aln.is_codon_aware():
            fr = aln.coding.frame
            keep = [c for c in range(aln.n_codons())
                    if not all(r[fr + 3 * c: fr + 3 * c + 3] == "---" for r in aln.rows)]
            rows = ["".join(r[fr + 3 * c: fr + 3 * c + 3] for c in keep) for r in aln.rows]
            new = Alignment(list(aln.ids), rows, aln.alphabet, coding=CodingSpec(0, aln.coding.table))
        else:
            keep = [c for c in range(aln.length) if any(r[c] not in "-." for r in aln.rows)]
            rows = ["".join(r[c] for c in keep) for r in aln.rows]
            new = Alignment(list(aln.ids), rows, aln.alphabet)
        removed = aln.length - new.length
        self._commit(new, f"removed {removed} gap-only column(s)")

    # ------------------------------------------------------------------ #
    def _commit(self, new_aln: Alignment, note: str):
        """Adopt a derived alignment, appending `note` to the provenance log so a
        hand-built (multi-program) alignment carries a reproducible recipe.

        This pushes onto the undo stack. Realigning discards every hand edit made
        since the last alignment, which is as destructive as any single edit and
        was previously the one destructive action with no way back.
        """
        hist = list(self.aln.meta.get("history", [])) if self.aln else []
        if note:
            hist.append(note)
        new_aln.meta["history"] = hist
        # Pins are indices into the alignment they were set on, so an edit that
        # changes the number of columns moves them. An operation that knows where
        # they went says so on the new alignment (anchored_realign does); one
        # that does not gets its pins carried only when the width is unchanged,
        # because a pin pointing at whatever slid into its place is worse than no
        # pin at all.
        if "anchors" not in new_aln.meta and self.aln is not None:
            if new_aln.length == self.aln.length:
                pins = list(self.aln.meta.get("anchors", ()))
                if pins:
                    new_aln.meta["anchors"] = pins
        # A residue mask names residues, not columns, and no edit or realignment
        # changes a sequence's residues, so it is carried whatever the width.
        # residue_mask() drops anything that no longer fits.
        if rel_mod.RESIDUE_MASK_KEY not in new_aln.meta and self.aln is not None:
            carried = rel_mod.residue_mask(rel_mod.with_residue_mask(
                new_aln, rel_mod.residue_mask(self.aln)))
            if carried:
                new_aln.meta[rel_mod.RESIDUE_MASK_KEY] = carried
        if self.aln is not None:
            self._undo_stack.append(self.aln)
            self._redo_stack.clear()
        self._set_alignment(new_aln)

    # ------------------------------------------------------------------ #
    # Truth mode
    # ------------------------------------------------------------------ #
    def _rescore_against_reference(self):
        """Kept as a name the window uses; the work is the document's."""
        self.doc.rescore()
        if self.doc.reference is None and self.track_combo.currentIndex() == _TRUTH_TRACK:
            self.track_combo.setCurrentIndex(0)
        self._truth_views_follow()

    # ----- showing the answer itself ---------------------------------------- #
    @property
    def showing_truth(self) -> bool:
        """Is the canvas displaying the reference rather than the document?"""
        return self._truth_view is not None

    def _toggle_true_alignment(self, on: bool):
        """Put the known-true alignment on screen, read-only, and take it off again.

        The document is not touched: nothing here changes the alignment, the
        undo stack or the dirty flag, and the canvas refuses edits for as long as
        the answer is up, because a keystroke aimed at the truth's columns would
        otherwise land on the working alignment's.
        """
        on = bool(on) and self.doc.truth_rows is not None
        if on == self.showing_truth:
            if not on:
                self._sync_truth_action(False)
            return
        if on:
            ref = Alignment(ids=list(self.aln.ids), rows=list(self.doc.truth_rows),
                            alphabet=self.aln.alphabet, coding=self.aln.coding,
                            meta={"source": "reference"})
            self._truth_view = {"level": self.canvas.level,
                                "colour": self.color_combo.currentIndex(),
                                "track": self.track_combo.currentIndex()}
            self.canvas.set_read_only(True)
            self.canvas.set_alignment(ref)
            self.canvas.set_column_scores(None)
            self.canvas.set_keep_mask(None)
            self.canvas.set_color_mode("residue")
            self.color_combo.blockSignals(True)
            self.color_combo.setCurrentIndex(_COLOUR_RESIDUE)
            self.color_combo.blockSignals(False)
            label = self._reference_label or "reference"
            self.sel_lbl.setText(f"showing the true alignment ({label}) — read-only; "
                                 "Ctrl+T returns to your alignment")
        else:
            state = self._truth_view or {}
            self._truth_view = None
            self.canvas.set_read_only(False)
            self.canvas.set_alignment(self.aln, level=state.get("level"))
            self.color_combo.blockSignals(True)
            self.color_combo.setCurrentIndex(state.get("colour", _COLOUR_RESIDUE))
            self.color_combo.blockSignals(False)
            self._on_colour()
            self.track_combo.setCurrentIndex(state.get("track", 0))
            self._report_accuracy()
        self._sync_truth_action(on)
        self._update_title()
        self._update_legend()
        self._refresh_inspector()
        self._refresh_controls()

    def _sync_truth_action(self, on: bool):
        self.show_truth_act.blockSignals(True)
        self.show_truth_act.setChecked(on)
        self.show_truth_act.blockSignals(False)

    def _truth_views_follow(self):
        """Keep the views that show the answer in step with whether there is one.

        The reference can go away underneath them — cleared by hand, or dropped
        because the document became a different set of sequences — and a stale
        answer on screen is worse than none.
        """
        if self.doc.truth_cells is not None:
            if self.color_combo.currentIndex() == _COLOUR_TRUTH:
                self._paint_truth_cells()
        else:
            if self.showing_truth:
                self._toggle_true_alignment(False)
            if self.color_combo.currentIndex() == _COLOUR_TRUTH:
                self.color_combo.setCurrentIndex(_COLOUR_RESIDUE)   # triggers _on_colour
        self._refresh_inspector()
        self._refresh_controls()

    def _generate_dataset(self):
        """Simulate a dataset, load the unaligned sequences, keep the truth."""
        if not self._confirm_discard("Generating a new dataset"):
            return
        dlg = SimulateDialog(self)
        if dlg.exec() != QDialog.Accepted:
            return
        v = dlg.values()
        from .. import simulate

        try:
            d = simulate.simulate(**v)
        except Exception as exc:                        # pragma: no cover - defensive
            QMessageBox.critical(self, "Could not generate a dataset", str(exc))
            return
        # Load the *unaligned* sequences: the whole exercise is to align them.
        # Right-padded, and flagged unaligned, the same way io.load_alignment
        # presents ragged input — an Alignment is rectangular by construction.
        records = list(d["seqs"])
        width = max(len(seq) for _, seq in records)
        unaligned = Alignment([n for n, _ in records],
                              [seq.ljust(width, "-") for _, seq in records],
                              Alphabet.DNA)
        unaligned.meta["unaligned"] = True
        unaligned.meta["history"] = [
            f"simulated: {v['taxa']} sequences, root {v['root_len']} nt, "
            f"divergence {v['bmax']}, indel rate {v['indel_rate']}, "
            f"rate alpha {v['rate_alpha']}, seed {v['seed']}"]
        self.doc.clear_reference()               # drop any previous reference first
        self._set_alignment(unaligned)
        self._annotations = []
        # …and hold the answer, scored straight away against the unaligned
        # sequences. Holding it without scoring left truth mode refusing to turn
        # on until something else changed the alignment, which reads as the
        # program having forgotten the dataset it just generated.
        self.doc.set_reference(
            Alignment(list(d["names"]), list(d["true_rows"]), Alphabet.DNA),
            "simulated truth")
        self._refresh_controls()
        self.sel_lbl.setText(
            f"{unaligned.n_seqs} simulated sequences loaded (unaligned). "
            "Press Align, then colour by Correct placement, or read the Column "
            "inspector, to see what the aligner got wrong and what it should "
            "have done.")

    def _load_reference(self):
        """Load a trusted alignment of the same sequences and turn on truth mode."""
        if self.aln is None:
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open reference alignment", "",
            "Alignments (*.fa *.fasta *.fas *.aln *.phy *.sto *.nex *.msf);;All files (*)")
        if not path:
            return
        try:
            ref = io.load_alignment(path, alphabet=self.aln.alphabet)
        except Exception as exc:
            title, text = _friendly_open_error(path, exc)
            QMessageBox.critical(self, title, text)
            return
        self._set_reference(ref, os.path.basename(path))

    def _set_reference(self, ref: Alignment, label: str):
        try:
            self.doc.set_reference(ref, label)
        except IdMismatch as exc:
            QMessageBox.critical(self, "Not a reference for this alignment", str(exc))
            return
        self.track_combo.setCurrentIndex(_TRUTH_TRACK)      # triggers _on_track
        self._report_accuracy()
        self._truth_views_follow()      # the answer is available now — offer it

    def _report_accuracy(self, label: Optional[str] = None):
        """Put the current accuracy in the status bar.

        Refreshed on every edit while the truth track is showing, so the numbers
        move as the alignment is worked on — which is what makes a hand edit
        legible as an improvement or a mistake rather than just a change.
        """
        if self.truth is None:
            return
        a = self.truth
        if label is None:
            label = self._reference_label
        auc = a.reliability_auc(self.reliability.col_combined) if self.reliability else None
        wrong = int(np.sum(np.nan_to_num(a.col_correct, nan=1.0) < 0.999))
        bits = [f"vs {label}" if label else "vs reference",
                f"SP={a.sp:.3f}", f"TC={a.tc:.3f}", f"precision={a.precision:.3f}",
                f"{wrong} columns wrong"]
        if auc is not None and auc == auc:
            bits.append(f"reliability AUC={auc:.3f}")
        self.sel_lbl.setText("  ".join(bits))

    def _clear_reference(self):
        if self.showing_truth:
            self._toggle_true_alignment(False)
        self.doc.clear_reference()
        if self.track_combo.currentIndex() == _TRUTH_TRACK:
            self.track_combo.setCurrentIndex(0)
        self._truth_views_follow()
        self._refresh_controls()

    # ------------------------------------------------------------------ #
    # Sessions, autosave and unsaved work
    # ------------------------------------------------------------------ #
    def _current_session(self, autosaved: bool = False) -> session_mod.Session:
        """Snapshot everything needed to resume this session. The document holds
        the data; the window contributes only the view state."""
        return self.doc.to_session(
            view={"track": self.track_combo.currentIndex(),
                  "threshold": self.thr.value(),
                  "level": self.canvas.level.value,
                  "colour": self.color_combo.currentIndex(),
                  "scheme": self.scheme_combo.currentText()},
            autosaved=autosaved)

    def _apply_session(self, sess: session_mod.Session, path: str = ""):
        """Restore a session. The view state is applied last, after the
        alignment, so the tracks it names have something to be computed from."""
        if sess.alignment is None:
            QMessageBox.warning(self, "Empty session", "That session holds no alignment.")
            return
        self.doc.restore_context(sess, path)
        self._set_alignment(sess.alignment, dirty=False)

        view = sess.view or {}
        if "threshold" in view:
            self.thr.blockSignals(True)
            self.thr.setValue(int(view["threshold"]))
            self.thr.blockSignals(False)
            self.thr_lbl.setText(f"{self.thr.value() / 100.0:.2f}")
        level = view.get("level")
        if level:
            for i in range(self.level_combo.count()):
                if getattr(self.level_combo.itemData(i), "value", None) == level:
                    self.level_combo.setCurrentIndex(i)
                    break
        scheme = view.get("scheme")
        if scheme and self.scheme_combo.findText(scheme) >= 0:
            self.scheme_combo.setCurrentText(scheme)
        track = view.get("track", 0)
        # An expensive track is not resumed automatically: restoring a session
        # should not silently start a minutes-long computation.
        if isinstance(track, int) and track in _CHEAP_TRACKS:
            self.track_combo.setCurrentIndex(track)
        self._set_dirty(False)

    # -- explicit session documents ---------------------------------------- #
    def _save_session(self, ask: bool = False):
        if self.aln is None:
            return False
        path = self._session_path
        if ask or not path:
            default = (session_mod.session_path_for(self._path) if self._path
                       else "session" + session_mod.SUFFIX)
            path, _ = QFileDialog.getSaveFileName(
                self, "Save session", default,
                f"CRAIC session (*{session_mod.SUFFIX});;All files (*)")
            if not path:
                return False
        try:
            self._current_session().save(path)
        except OSError as exc:
            QMessageBox.critical(self, "Could not save the session", str(exc))
            return False
        self._session_path = path
        self._set_dirty(False)
        self._clear_autosave()
        self.sel_lbl.setText(f"session saved to {os.path.basename(path)}")
        return True

    def _open_session(self):
        if not self._confirm_discard("Opening another session"):
            return
        path, _ = QFileDialog.getOpenFileName(
            self, "Open session", "",
            f"CRAIC session (*{session_mod.SUFFIX} *.json);;All files (*)")
        if path:
            self.open_session_path(path)

    def open_session_path(self, path: str) -> bool:
        try:
            sess = session_mod.Session.load(path)
        except session_mod.SessionError as exc:
            QMessageBox.critical(self, "Could not open the session", str(exc))
            return False
        self._apply_session(sess, path)
        return True

    # -- crash recovery ----------------------------------------------------- #
    def _autosave_path(self) -> str:
        """Where the recovery file lives: the application's own state directory,
        never beside the user's data. Autosaving next to an input file would
        mean that opening someone's reference alignment and looking at it left
        files behind in their directory."""
        from PySide6.QtCore import QStandardPaths

        base = QStandardPaths.writableLocation(QStandardPaths.AppDataLocation) or ""
        if not base:                              # pragma: no cover - defensive
            base = os.path.join(os.path.expanduser("~"), ".craic")
        return os.path.join(base, "craic", "recovery" + session_mod.SUFFIX)

    def _touch_autosave(self):
        """Schedule an autosave shortly after a change, coalescing bursts of
        edits into one write."""
        if self.aln is None:
            return
        self._autosave_timer.start()

    def _write_autosave(self):
        if self.aln is None or not self._dirty:
            return
        try:
            self._current_session(autosaved=True).save(self._autosave_path())
        except OSError:
            pass          # best effort: never interrupt the user to report this

    def _clear_autosave(self):
        """Remove the recovery file. Called on a clean exit and on an explicit
        save, so that the presence of one always means something went wrong."""
        try:
            os.remove(self._autosave_path())
        except OSError:
            pass

    def maybe_recover(self) -> bool:
        """Offer to restore an interrupted session. Returns True if one was."""
        path = self._autosave_path()
        if not os.path.exists(path):
            return False
        try:
            sess = session_mod.Session.load(path)
        except session_mod.SessionError:
            self._clear_autosave()
            return False
        when = time.strftime("%H:%M on %d %b", time.localtime(sess.saved_at))
        answer = QMessageBox.question(
            self, "Recover unsaved work",
            f"CRAIC closed without saving.\n\n{sess.describe()}\n\n"
            f"Last change at {when}. Recover it?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        if answer != QMessageBox.Yes:
            self._clear_autosave()
            return False
        self._apply_session(sess)
        self._set_dirty(True)          # recovered, but still not saved anywhere
        return True

    # -- replacing the document --------------------------------------------- #
    def _confirm_discard(self, what: str) -> bool:
        """Ask before throwing unsaved work away. Returns False to abandon.

        Every way of losing a session goes through here: quitting, opening
        another file, opening another session, generating a dataset. They are the
        same question, so they are the same prompt — the previous version had two
        copies with different wording, and adding a third place that needed it is
        how an aligned dataset got discarded silently.
        """
        if not self.doc.dirty or self.doc.is_empty:
            return True
        answer = QMessageBox.warning(
            self, "Unsaved changes",
            "This alignment has changes that have not been saved.\n\n"
            f"{what} will discard them.\n\n"
            "Saving the session keeps everything — the alignment, the reference, "
            "annotations and the provenance log. Saving the alignment alone "
            "writes a standard format and keeps only the residues.",
            QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
            QMessageBox.Save)
        if answer == QMessageBox.Cancel:
            return False
        if answer == QMessageBox.Save:
            return self._save_session()
        return True

    # -- quitting ----------------------------------------------------------- #
    def closeEvent(self, event):
        """Never discard a curation session silently."""
        if self._confirm_discard("Quitting"):
            self._stop_job()
            self._clear_autosave()
            event.accept()
        else:
            event.ignore()

    # ------------------------------------------------------------------ #
    # The document, under the names the rest of the window (and the tests)
    # already use. These are deliberately thin: the state lives in self.doc.
    # ------------------------------------------------------------------ #
    @property
    def aln(self):
        return self.doc.alignment

    @property
    def reference(self):
        return self.doc.reference

    @reference.setter
    def reference(self, value):
        self.doc._reference = value
        if value is None:
            self.doc.reference_label = ""

    @property
    def truth(self):
        return self.doc.truth

    @truth.setter
    def truth(self, value):
        self.doc.truth = value

    @property
    def _reference_label(self):
        return self.doc.reference_label

    @_reference_label.setter
    def _reference_label(self, value):
        self.doc.reference_label = value

    @property
    def _dirty(self):
        return self.doc.dirty

    @property
    def _groups(self):
        return self.doc.groups

    @_groups.setter
    def _groups(self, value):
        self.doc.groups = value

    @property
    def _annotations(self):
        return self.doc.annotations

    @_annotations.setter
    def _annotations(self, value):
        self.doc.annotations = value

    @property
    def _path(self):
        return self.doc.path

    @_path.setter
    def _path(self, value):
        self.doc.path = value

    @property
    def _session_path(self):
        return self.doc.session_path

    @_session_path.setter
    def _session_path(self, value):
        self.doc.session_path = value

    @property
    def _doc_generation(self):
        return self.doc.generation

    def _set_alignment(self, aln: Alignment, dirty: bool = True):
        # ``dirty`` is False only when the document is being *replaced* rather
        # than changed — opening a file, recovering a session — since those are
        # not unsaved work.
        # A cheap overlay is kept across the change and recomputed below, so that
        # editing the alignment updates the track rather than clearing it. This
        # matters most for truth mode: watching correctness change as you move a
        # block is the whole point of having the answer.
        previous_track = self.track_combo.currentIndex()
        previous_colour = self.color_combo.currentIndex()
        self.doc.set_alignment(aln, dirty=dirty)      # also rescores against the reference
        self._stop_job()                   # an analysis of the old alignment is no use now
        self.reliability = None
        self.agreement = None
        self._reliability_waiters = []
        self._perturbation_skipped = False
        self._consistency_secs = None
        self._warned_frame_breaks = False     # warn again for this new alignment
        self.canvas.set_alignment(aln)
        self.canvas.set_annotations(self._annotations)
        self.canvas.set_pins(aln.meta.get("anchors", ()))
        self.canvas.set_residue_mask(rel_mod.masked_cells(aln, rel_mod.residue_mask(aln)))
        self._sync_engine_combo()
        self.canvas.set_column_scores(None)
        self.canvas.set_keep_mask(None)
        self.posterior.set_alignment(aln)
        self.sandbox.set_context(aln, self._engines)
        # rebuild level combo
        self.level_combo.blockSignals(True)
        self.level_combo.clear()
        for lv in aln.available_levels():
            self.level_combo.addItem(_LEVEL_LABEL[lv], lv)
        idx = max(0, [self.level_combo.itemData(i) for i in range(self.level_combo.count())].index(self.canvas.level))
        self.level_combo.setCurrentIndex(idx)
        self.level_combo.blockSignals(False)
        self.coding_chk.blockSignals(True)
        self.coding_chk.setChecked(aln.coding is not None)
        self.coding_chk.setEnabled(aln.alphabet.is_nucleotide)
        self.coding_chk.blockSignals(False)
        self._sync_code_combo()
        keep_track = previous_track if previous_track in _CHEAP_TRACKS else 0
        if keep_track == _TRUTH_TRACK and self.truth is None:
            keep_track = 0                     # the reference no longer applies
        self.track_combo.blockSignals(True)
        self.track_combo.setCurrentIndex(keep_track)
        self.track_combo.blockSignals(False)
        self.outlier_chk.blockSignals(True)
        self.outlier_chk.setChecked(False)
        self.outlier_chk.blockSignals(False)
        self.canvas.set_outlier_rows(set())
        # Per-residue truth colouring survives the change for the same reason the
        # truth track does: seeing which residues an edit put right is the point.
        # Everything else goes back to plain residue colours, since the analysis
        # behind it was computed for the alignment that has just been replaced.
        keep_colour = (_COLOUR_TRUTH if previous_colour == _COLOUR_TRUTH
                       and self.doc.truth_cells is not None else _COLOUR_RESIDUE)
        self.color_combo.blockSignals(True)
        self.color_combo.setCurrentIndex(keep_colour)
        self.color_combo.blockSignals(False)
        if keep_colour == _COLOUR_TRUTH:
            self._paint_truth_cells()
        else:
            self.canvas.set_cell_scores(None)
            self.canvas.set_color_mode("residue")
        self.find_edit.clear()
        if self._inspect_col is not None and self._inspect_col >= aln.length:
            self._inspect_col = None       # the alignment shrank out from under it
        self._refresh_inspector()
        self._update_tool_level()
        self._update_legend()
        self._refresh_controls()
        if self.track_combo.currentIndex() != 0:
            self._on_track()                   # recompute the kept overlay
        self._set_dirty(dirty)

    def _set_dirty(self, dirty: bool = True):
        self.doc.set_dirty(dirty)
        self._update_title()
        if dirty:
            self._touch_autosave()

    def _update_title(self):
        name = os.path.basename(self._path) if self._path else ""
        if not name and self._session_path:
            name = os.path.basename(self._session_path)
        base = "CRAIC — Conserved-Region Alignment by Iterative Convergence"
        if name:
            base = f"{name} — CRAIC"
        if self.showing_truth:
            base = f"{base} — true alignment (read-only)"
        self.setWindowTitle(("• " if self._dirty else "") + base)
        # macOS draws its own close-button dot from this; harmless elsewhere.
        self.setWindowModified(False)

    #: The engine selected at start-up. The built-in one, deliberately: it is the
    #: only engine guaranteed to be present, so it is the only one that can be a
    #: default on every machine, and it is the engine whose posteriors the
    #: reliability overlay, the posterior explorer, the sandbox and the
    #: perturbation ensemble all read — so starting there makes the workbench
    #: coherent before the user chooses anything. On BAliBASE core blocks it is
    #: level with MAFFT (0.835 vs 0.830 SP, n=188, Wilcoxon p=0.51) and ahead of
    #: Clustal Omega, though behind MUSCLE and ProbCons; the Engine menu is one
    #: click away for work where that matters.
    _DEFAULT_ENGINE = "builtin"

    def _default_engine_index(self) -> int:
        """Index of the start-up engine, or 0 if it somehow is not registered."""
        for i, e in enumerate(self._engines):
            if e.key == self._DEFAULT_ENGINE:
                return i
        return 0

    def _engine_alphabet(self) -> Alphabet:
        """The alphabet the engine will actually be handed.

        With *align as protein* ticked, coding nucleotides are translated before
        the engine sees them, so a protein-only engine is perfectly usable — it
        is handed amino acids. Asking about the document's own alphabet instead
        would rule out exactly the case the checkbox exists to enable.
        """
        if (self.aln is not None and self.aln.alphabet.is_nucleotide
                and self.codon_chk.isChecked()):
            return Alphabet.PROTEIN
        return self.aln.alphabet if self.aln else Alphabet.PROTEIN

    def _sync_engine_combo(self):
        """Disable engines that cannot align the data as it will reach them.

        A protein-only engine in the list is a trap unless the list says so: the
        entry stays visible, so the user can see the tool exists and why it is
        unavailable, rather than wondering where it went.
        """
        if not self.aln:
            return
        alphabet = self._engine_alphabet()
        model = self.engine_combo.model()
        for i in range(self.engine_combo.count()):
            eng = self.engine_combo.itemData(i)
            ok = eng is None or eng.supports(alphabet)
            item = model.item(i) if hasattr(model, "item") else None
            if item is not None:
                item.setEnabled(ok)
            self.engine_combo.setItemData(
                i, "" if ok else (f"{eng.label} aligns protein sequences only — "
                                  "tick “align as protein” to translate first"),
                Qt.ToolTipRole)
        cur = self.engine_combo.currentData()
        if cur is not None and not cur.supports(alphabet):
            self.engine_combo.setCurrentIndex(self._default_engine_index())

    def _refresh_controls(self):
        has = self.aln is not None
        showing = self.showing_truth
        self.align_btn.setEnabled(has and not showing)
        for a in self._doc_actions:           # Save/Copy/Export/Sort/… need an alignment
            a.setEnabled(has)
        for a in self._edit_actions:          # nothing that rewrites the alignment
            a.setEnabled(has and not showing)
        self.track_combo.setEnabled(has and not showing)
        self.color_combo.setEnabled(has and not showing)
        self.show_truth_act.setEnabled(self.doc.truth_rows is not None)
        if has and showing:
            return                            # the status line says what is on screen
        if has:
            self.sel_lbl.setText(
                f"{self.aln.n_seqs} seqs × {self.aln.length} cols ({self.aln.alphabet.value})")

    def _update_tool_level(self):
        """Keep the sandbox and posterior explorer working at the viewed level."""
        self._update_scheme_enabled()
        if not self.aln:
            return
        level = self.canvas.level
        coding = self.aln.coding is not None
        table = self.aln.coding.table if self.aln.coding else 1
        aa_mode = coding and level in (Level.CODON, Level.AA)
        self.sandbox.set_view(aa_mode, table, level)
        self.posterior.set_view(aa_mode)

    def _update_scheme_enabled(self):
        """Grey out the residue-colour scheme when it would do nothing.

        All the schemes are amino-acid schemes, so in a nucleotide view they have
        no effect; leaving the control live there means changing it appears to
        fail.
        """
        applies = self.aln is not None and colors.scheme_applies(
            self.canvas.level, self.aln.alphabet)
        self.scheme_combo.setEnabled(applies)
        self.scheme_lbl.setEnabled(applies)
        self.scheme_combo.setToolTip(
            "Residue colour scheme."
            if applies else
            "The colour schemes are amino-acid schemes. Nucleotides are coloured "
            "by base, so this has no effect in the nucleotide view — switch View "
            "to Codon or Amino acid (or open a protein alignment) to use it.")

    def _on_level(self, *_):
        lv = self.level_combo.currentData()
        if lv is not None:
            self.canvas.set_level(lv)
            self._update_tool_level()
            if lv in (Level.CODON, Level.AA):
                self._maybe_warn_frame_breaks()

    def _maybe_warn_frame_breaks(self):
        """When showing a coding alignment as codon/protein, flag codons that
        straddle a gap (they show as X) and offer the codon-aware fix. Once per
        alignment so it isn't naggy."""
        if not self.aln or self.aln.coding is None or self._warned_frame_breaks:
            return
        breaks = self.aln.frame_break_columns()
        if not breaks:
            return
        self._warned_frame_breaks = True
        n = len(breaks)
        self.sel_lbl.setText(f"{n} codon(s) straddle a gap and show as X")
        QMessageBox.warning(
            self, "Reading frame not codon-aligned",
            f"{n} codon column(s) in this alignment straddle a gap, so they can't be "
            "translated and show as X in the codon / amino-acid view.\n\n"
            "This happens when the gaps weren't placed on codon boundaries. To get a "
            "clean protein view, either re-align with “align as protein” ticked, or "
            "select an affected region and apply a (codon-aware) alternative from the "
            "Realignment sandbox.")

    def _sync_code_combo(self):
        """Show the document's genetic code, and offer it only when it applies."""
        table = self.aln.coding.table if (self.aln and self.aln.coding) else 1
        self.code_combo.blockSignals(True)
        idx = self.code_combo.findData(table)
        self.code_combo.setCurrentIndex(idx if idx >= 0 else 0)
        self.code_combo.blockSignals(False)
        self.code_combo.setEnabled(
            bool(self.aln) and self.aln.alphabet.is_nucleotide)

    def _current_code(self) -> int:
        return self.code_combo.currentData() or 1

    def _on_code(self, *_):
        """Re-read the alignment under a different genetic code."""
        if not self.aln or not self.aln.coding:
            return
        self.aln.coding = CodingSpec(frame=self.aln.coding.frame,
                                     table=self._current_code())
        self.canvas.set_alignment(self.aln, level=self.canvas.level)
        self._update_tool_level()
        self._warn_internal_stops()

    def _warn_internal_stops(self):
        """Tell the user when a coding alignment stops in mid-sequence.

        Usually the wrong genetic code, the wrong reading frame, or a
        pseudogene. CRAIC cannot tell which, but saying nothing is the worst
        option: under the standard code a mitochondrial or Mollicute gene fills
        with stop codons, and a stop carries no homology signal at all, so
        “align as protein” quietly loses exactly the conserved tryptophans an
        aligner most relies on.
        """
        if not self.aln or not self.aln.coding:
            return
        stops = self.aln.internal_stops()
        if not stops:
            return
        rows = sorted({r for r, _ in stops})
        table = self.aln.coding.table
        names = ", ".join(self.aln.ids[r] for r in rows[:4])
        if len(rows) > 4:
            names += f" and {len(rows) - 4} more"
        extra = ("\n\nMitochondrial and Mollicute (Mycoplasma, Spiroplasma) "
                 "sequences read TGA as tryptophan. Under the standard code every "
                 "one of those becomes a stop, which contributes no homology "
                 "signal — try genetic code 2 or 4."
                 if table == 1 else
                 "\n\nCheck the reading frame, or whether these are pseudogenes.")
        self.sel_lbl.setText(f"{len(stops)} internal stop codons "
                             f"({code_name(table)} code)")
        QMessageBox.warning(
            self, "Internal stop codons",
            f"{len(stops)} stop codons fall inside {len(rows)} of "
            f"{self.aln.n_seqs} sequences ({names}), reading the "
            f"{code_name(table)} code.{extra}")

    def _on_coding(self, *_):
        if not self.aln:
            return
        turned_on = self.coding_chk.isChecked() and self.aln.alphabet.is_nucleotide
        if turned_on:
            self.aln.coding = CodingSpec(frame=0, table=self._current_code())
        else:
            self.aln.coding = None
        # rebuild levels, keep scores
        self.canvas.set_alignment(self.aln)
        self.level_combo.blockSignals(True)
        self.level_combo.clear()
        for lv in self.aln.available_levels():
            self.level_combo.addItem(_LEVEL_LABEL[lv], lv)
        self.level_combo.blockSignals(False)
        self._update_tool_level()
        self._sync_code_combo()
        self._on_track()
        if turned_on:
            self._warn_internal_stops()
        if self.canvas.level in (Level.CODON, Level.AA):
            self._maybe_warn_frame_breaks()

    # ------------------------------------------------------------------ #
    def _source_records(self):
        return [(i, r.replace("-", "")) for i, r in zip(self.aln.ids, self.aln.rows)]

    def _params_for(self, engine) -> dict:
        defaults = {p.key: p.default for p in engine.parameters()}
        defaults.update(self._engine_params.get(engine.key, {}))
        return defaults

    def _edit_engine_params(self):
        engine = self.engine_combo.currentData()
        if engine is None:
            return
        if not engine.parameters():
            QMessageBox.information(self, engine.label,
                                    "This aligner has no tunable parameters in CRAIC.")
            return
        dlg = EngineParamsDialog(engine, self._params_for(engine), self,
                                 alphabet=self._engine_alphabet() if self.aln else None)
        if dlg.exec():
            self._engine_params[engine.key] = dlg.values()
            self.statusBar().showMessage(f"{engine.label} parameters updated", 3000)

    def _align(self):
        if not self.aln:
            return
        engine = self.engine_combo.currentData()
        params = self._params_for(engine)
        self._want_aa = self.codon_chk.isChecked() and self.aln.alphabet.is_nucleotide
        # Refuse an engine that cannot take this data *here*, where we can explain
        # it and offer the way round, rather than letting it raise inside a worker
        # thread and surface as a traceback.
        if not engine.supports(self._engine_alphabet()):
            extra = ("\n\nThese are coding nucleotides, so ticking “align as protein” "
                     "will translate them, align the amino acids and thread the "
                     "nucleotides back through the result."
                     if self.aln.coding is not None else
                     "\n\nChoose a different engine, or tick “treat as coding” and "
                     "then “align as protein” if these are protein-coding sequences.")
            QMessageBox.information(
                self, f"{engine.label} cannot align this data",
                f"{engine.label} aligns protein sequences only, and this alignment "
                f"is {self.aln.alphabet.value}." + extra)
            return
        if self._want_aa:
            engine = CodonAware(engine, table=self._current_code())
        # the History lists only settings that mean something for this data
        handed = self._engine_alphabet()
        shown = _fmt_params({p.key: params[p.key] for p in engine.parameters()
                             if p.applies_to(handed) and p.key in params})
        self._align_note = engine.label + (f" ({shown})" if shown else "")
        records = self._source_records()
        alphabet = self.aln.alphabet
        self.align_btn.setEnabled(False)
        self._align_cancelled = False
        self._align_base = (f"{engine.label}\n"
                            f"{len(records):,} sequences × {self.aln.length:,} columns")

        dlg = QProgressDialog(self._align_base + "…", "Cancel", 0, 0, self)
        dlg.setWindowTitle("Aligning")
        dlg.setWindowModality(Qt.WindowModal)
        dlg.setMinimumDuration(0)
        dlg.setAutoClose(False)
        dlg.setAutoReset(False)
        dlg.canceled.connect(self._cancel_align)
        self._align_progress = dlg

        def work(report, cancelled):
            return engine.align(records, alphabet, progress=report, cancelled=cancelled, **params)

        self._align_token = run_cancellable(work, self._on_aligned, self._on_align_error,
                                            self._on_align_progress)
        dlg.show()

    def _on_align_progress(self, done: int, total: int, label: str):
        dlg = self._align_progress
        if dlg is None:
            return
        if total > 0 and dlg.maximum() != total:
            dlg.setRange(0, total)
        if not self._align_cancelled:
            base = getattr(self, "_align_base", "Aligning")
            dlg.setLabelText(f"{base}\nmerging profiles {done} of {total}…"
                             if total else base + "…")
        dlg.setValue(done)        # modal dialog re-enters the loop here — do it last

    def _cancel_align(self):
        self._align_cancelled = True
        if self._align_token is not None:
            self._align_token.cancel()
        if self._align_progress is not None:
            self._align_progress.setLabelText("Cancelling…")

    def _close_align_progress(self):
        if self._align_progress is not None:
            try:
                self._align_progress.canceled.disconnect(self._cancel_align)
            except (RuntimeError, TypeError):
                pass
            self._align_progress.close()
            self._align_progress = None
        self._align_token = None

    def _on_aligned(self, aln):
        self._close_align_progress()
        self._idle()
        if self._align_cancelled:
            self._align_cancelled = False
            return
        self._commit(aln, f"aligned · {self._align_note}")
        # "Align as protein" should land you on the protein view.
        if getattr(self, "_want_aa", False) and Level.AA in aln.available_levels():
            idx = next((i for i in range(self.level_combo.count())
                        if self.level_combo.itemData(i) == Level.AA), -1)
            if idx >= 0:
                self.level_combo.setCurrentIndex(idx)   # triggers _on_level → AA view
        self._warn_internal_stops()

    def _on_align_error(self, msg):
        self._close_align_progress()
        self._idle()
        if self._align_cancelled:
            self._align_cancelled = False
            return
        # The worker hands us a full traceback. A traceback is what a developer
        # needs and the last thing a user needs, so lead with the exception's own
        # message and keep the rest behind "Show Details".
        lines = [ln for ln in (msg or "").strip().splitlines() if ln.strip()]
        summary = lines[-1] if lines else "The aligner failed."
        if ": " in summary:
            summary = summary.split(": ", 1)[1]
        box = QMessageBox(QMessageBox.Warning, "Alignment failed", summary, QMessageBox.Ok, self)
        box.setInformativeText("The alignment on screen has not been changed.")
        box.setDetailedText(msg[-4000:])
        box.exec()

    # ------------------------------------------------------------------ #
    # Long analyses: one at a time, off the GUI thread, behind a progress
    # dialog that shows the time left and can cancel.
    # ------------------------------------------------------------------ #
    def _start_job(self, kind, what, work, on_done, on_cancel=None, on_error=None):
        """Run ``work(report, cancelled)`` in the background.

        ``on_done(result, seconds)`` runs when it finishes, ``on_cancel()`` if the
        user cancels it and ``on_error(message)`` if it fails, each only while
        the same alignment is on screen: a result for a replaced alignment is
        dropped.
        """
        aln = self.aln
        text = f"{what}\n{aln.n_seqs:,} sequences × {aln.length:,} columns"
        if accel.backend() != "rust":
            text += "\n\n" + _NO_CORE
        dlg = JobProgress("CRAIC", text, self)
        job = {"kind": kind, "dialog": dlg, "cancelled": False, "token": None}
        gen, t0 = self._doc_generation, time.monotonic()

        def finish(callback, *args):
            _close_job_dialog(job)
            if self._job is job:
                self._job = None
                self._idle_track()
            if gen != self._doc_generation:
                return
            if job["cancelled"]:
                callback, args = on_cancel, ()
            if callback is not None:
                callback(*args)

        def progress(done, total, label):
            if self._job is job:
                dlg.update_progress(done, total, label)

        dlg.canceled.connect(lambda: self._cancel_job(job))
        job["token"] = run_cancellable(
            work, lambda r: finish(on_done, r, time.monotonic() - t0),
            lambda msg: finish(on_error or self._on_track_error, msg), progress)
        self._job = job
        self._busy_track(what + "…")

    def _cancel_job(self, job):
        job["cancelled"] = True
        job["token"].cancel()
        self.statusBar().showMessage("Cancelling…")

    def _stop_job(self):
        """Cancel the running analysis, if any, without waiting for it."""
        job, self._job = self._job, None
        if job is not None:
            self._cancel_job(job)
            _close_job_dialog(job)
            self._idle_track()

    def _busy_elsewhere(self, kind) -> bool:
        """True, with a word to the user, if a different analysis is running."""
        if self._job is not None and self._job["kind"] != kind:
            self.statusBar().showMessage(
                "Another analysis is running — wait for it to finish, or cancel it.", 6000)
            return True
        return False

    def _reset_overlay(self, kind):
        """Put the track (and, for reliability, the colouring) back to plain when
        the analysis it was waiting for is cancelled."""
        track = _TRACK_DEFS[self.track_combo.currentIndex()]
        if track.kind == kind:
            self.track_combo.blockSignals(True)
            self.track_combo.setCurrentIndex(0)
            self.track_combo.blockSignals(False)
            self.canvas.set_column_scores(None)
        if kind == "reliability" and self.color_combo.currentIndex() == _COLOUR_CONFIDENCE:
            self.color_combo.blockSignals(True)
            self.color_combo.setCurrentIndex(_COLOUR_RESIDUE)
            self.color_combo.blockSignals(False)
            self.canvas.set_color_mode("residue")
        self._update_legend()

    #: Ask before starting a perturbation ensemble predicted to take longer (s).
    _ASK_ABOVE_S = 120

    def _ensure_reliability(self, then):
        """Run ``then`` once the reliability report is ready.

        Computed once per alignment, in the background: the consistency score
        first, which is quick, then the perturbation ensemble, which re-aligns
        the sequences sixteen times and is nearly all of the cost. Every
        sequence pair costs about the same in both, so the consistency stage's
        time predicts the ensemble's; if that is long the user is asked first,
        and can settle for the consistency score alone.
        """
        rep = self.reliability
        if rep is not None and (rep.n_replicates_ok or self._perturbation_skipped):
            self.canvas.set_cell_scores(rep.cell_combined)
            then()
            return
        if self._busy_elsewhere("reliability"):
            return
        self._reliability_waiters.append(then)
        if self._job is not None:
            return                                   # already on its way
        if rep is None:
            self._run_consistency()
        else:
            self._offer_perturbation()

    def _run_consistency(self):
        aln = self.aln

        def work(report, cancelled):
            return rel_mod.consistency(aln, progress=report, cancelled=cancelled)

        def done(cons, secs):
            self._consistency_secs = secs
            self.reliability = rel_mod.combine(aln, cons)
            self._offer_perturbation()

        self._start_job("reliability", "Reliability: consistency score", work, done,
                        on_cancel=self._drop_waiters, on_error=self._reliability_failed)

    def _offer_perturbation(self):
        aln, rep = self.aln, self.reliability
        n_c, n_p = rel_mod.work(aln)
        # a pair costs the same in both stages; merging adds a little on top
        est = 1.2 * (self._consistency_secs or 0.0) / max(1, n_c) * n_p
        if est > self._ASK_ABOVE_S and not self._confirm_perturbation(est):
            self._skip_perturbation()
            return
        cons = (rep.col_consistency, rep.cell_consistency)

        def work(report, cancelled):
            stats: dict = {}
            pert = rel_mod.perturbation(aln, stats=stats, progress=report, cancelled=cancelled)
            return pert, stats

        def done(result, _secs):
            pert, stats = result
            self.reliability = rel_mod.combine(aln, cons, pert, stats)
            self._release_waiters()

        self._start_job("reliability", "Reliability: perturbation score (16 re-alignments)",
                        work, done, on_cancel=self._skip_perturbation,
                        on_error=self._reliability_failed)

    def _confirm_perturbation(self, seconds: float) -> bool:
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Question)
        box.setWindowTitle("Perturbation score")
        box.setText("The consistency score is ready. The perturbation score will take "
                    f"about {duration(seconds)} more on this computer.")
        box.setInformativeText(
            f"It re-aligns all {self.aln.n_seqs:,} sequences sixteen times. You can cancel "
            "it part-way and keep the consistency score. Without it, Reliability shows "
            "the consistency score alone."
            + ("\n\n" + _NO_CORE if accel.backend() != "rust" else ""))
        run = box.addButton("Run it", QMessageBox.AcceptRole)
        box.addButton("Consistency only", QMessageBox.RejectRole)
        box.setDefaultButton(run)
        box.exec()
        return box.clickedButton() is run

    def _skip_perturbation(self):
        self._perturbation_skipped = True
        self.statusBar().showMessage(
            "Reliability is the consistency score alone: the perturbation score was not run.",
            8000)
        self._release_waiters()

    def _release_waiters(self):
        waiters, self._reliability_waiters = self._reliability_waiters, []
        if self.reliability is not None:
            self.canvas.set_cell_scores(self.reliability.cell_combined)
        for then in waiters:
            then()

    def _drop_waiters(self):
        self._reliability_waiters = []
        self._reset_overlay("reliability")

    def _reliability_failed(self, msg):
        self._reliability_waiters = []
        self._on_track_error(msg)

    def _on_colour(self, *_):
        idx = self.color_combo.currentIndex()
        if idx == _COLOUR_CONFIDENCE:
            def show():
                self.canvas.set_color_mode("confidence")
                self._update_legend()
            self._ensure_reliability(show)
        elif idx == _COLOUR_TRUTH:
            if not self._paint_truth_cells():
                self.color_combo.blockSignals(True)
                self.color_combo.setCurrentIndex(_COLOUR_RESIDUE)
                self.color_combo.blockSignals(False)
                self.canvas.set_color_mode("residue")
                QMessageBox.information(self, "No known answer", self._no_truth_message())
            self._update_legend()
        else:
            self.canvas.set_color_mode("residue")
            self._update_legend()

    def _no_truth_message(self) -> str:
        """Why truth mode is unavailable, in terms of what to do about it.

        The two cases read identically on screen but need opposite actions, and
        telling a user to load a reference when the document is already holding
        one reads as the program having forgotten it.
        """
        if self.doc.reference is not None:
            label = self._reference_label or "a reference"
            return (f"The answer is held ({label}), but nothing has been scored "
                    "against it yet.\n\nPress Align — correctness is measured "
                    "against the reference every time the alignment changes.")
        return ("Colouring by correct placement needs an alignment whose true "
                "answer is known.\n\nEither generate a dataset "
                "(Teach ▸ Generate dataset with known answer…), which keeps its "
                "own truth, or load a trusted alignment of the same sequences "
                "(Teach ▸ Load reference alignment…).")

    def _paint_truth_cells(self) -> bool:
        """Push per-residue correctness onto the canvas. False if there is no truth.

        Called again after every edit while this mode is on, so that dragging a
        block recolours the residues it moved — the whole point of holding the
        answer is watching it respond.
        """
        cells = self.doc.truth_cells
        if cells is None:
            return False
        self.canvas.set_cell_scores(cells)
        self.canvas.set_color_mode("truth")
        return True

    def _update_legend(self):
        """Show the reliability colour key exactly when the ramp is on screen."""
        mode = self.canvas.color_mode
        painted = mode in ("confidence", "truth") or self.canvas.nt_scores is not None
        self.legend.set_meaning(mode)
        self.legend.setVisible(bool(self.aln is not None and painted))

    def _on_residue_clicked(self, seq: int, nt_col: int):
        if not self.aln:
            return
        row = self.aln.rows[seq]
        if nt_col >= len(row) or row[nt_col] == "-":
            self.canvas.set_probe(None)
            return
        res_index = sum(1 for ch in row[:nt_col] if ch != "-")
        self._last_probe = (seq, res_index)        # remembered for the homology-card export
        aln = self.aln
        model = post_mod.model_for(aln)
        self.statusBar().showMessage(f"Probing residue {res_index + 1} of {aln.ids[seq]}…")
        run_async(lambda: post_mod.probe_residue(aln, seq, res_index, model),
                  lambda probs: self._on_probe(probs, seq, nt_col),
                  self._on_error)

    def _on_probe(self, probs, seq, nt_col):
        self.statusBar().clearMessage()
        self.canvas.set_probe(probs, seq, nt_col)
        self.sel_lbl.setText("probe: brighter columns = where this residue could also align")

    # ----- figures, reports, hotspots --------------------------------------- #
    def _region_nt(self):
        nt = self.canvas.selection_nt()
        return nt if nt else (0, self.aln.length)

    def _render_fig(self, v: dict, path: str):
        """Render a figure from an explicit parameter dict (see FigureExportDialog)."""
        aln = self.aln
        level = self.canvas.level
        rel = self.reliability.col_combined if self.reliability else None
        theme = v.get("theme", "Light")
        sel = self.canvas.selection_nt()
        if v.get("selection") and sel:
            c0, c1 = sel
        else:
            c0, c1 = 0, aln.length
        kind = v["kind"]
        if kind == "alignment":
            opt = figures.FigureOptions(level=level, scheme=self.canvas.scheme,
                                        reliability=rel, title="CRAIC — alignment", theme=theme)
            figures.export_alignment(aln, path, opt=opt)
        elif kind == "report":
            figures.confidence_report(aln, self.reliability.col_combined, path, theme=theme)
        elif kind == "logo":
            figures.uncertainty_logo(aln, c0, c1, path, reliability=rel, level=level, theme=theme)
        elif kind == "arcs":
            figures.homology_arcs(aln, v["seq_i"], v["seq_j"], c0, c1, path, theme=theme)
        elif kind == "cooc":
            engines = self._engines if len(self._engines) > 1 else builtin_variants()
            co, ci, cj, n = ensemble_mod.cooccurrence(
                aln, v["seq_i"], v["seq_j"], c0, c1, engines,
                params_by_key={e.key: self._params_for(e) for e in engines})
            figures.matrix_heatmap(co, list(ci), list(cj), path, theme=theme,
                                   title="Ensemble co-occurrence",
                                   subtitle=f"{aln.ids[v['seq_i']]} vs {aln.ids[v['seq_j']]} · "
                                            f"{n} alternative alignments")
        elif kind == "card":
            figures.homology_card(aln, v["seq_i"], v["residue"], path, theme=theme)
        elif kind == "html":
            figures.export_html(aln, path, reliability=rel, level=level, scheme=self.canvas.scheme)

    def _export_fig(self, *_):
        if not self.aln:
            return
        probe = getattr(self, "_last_probe", None)
        dlg = FigureExportDialog(
            self.aln, self.posterior.seq_i.currentIndex(), self.posterior.seq_j.currentIndex(),
            probe, self.canvas.selection_nt() is not None, self)
        if not dlg.exec():
            return
        v = dlg.values()
        ext = {"pdf": ".pdf", "svg": ".svg", "png": ".png", "html": ".html"}[v["fmt"]]
        filt = {"pdf": "PDF (*.pdf)", "svg": "SVG (*.svg)", "png": "PNG image (*.png)",
                "html": "HTML (*.html)"}[v["fmt"]]
        path, _ = QFileDialog.getSaveFileName(self, "Export figure", v["kind"] + ext, filt)
        if not path:
            return
        if v["kind"] == "report":          # needs the reliability report: wait for it
            self._ensure_reliability(lambda: self._write_fig(v, path))
        else:
            self._write_fig(v, path)

    def _write_fig(self, v: dict, path: str):
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            self._render_fig(v, path)
            self.statusBar().showMessage(f"Saved figure → {path}", 5000)
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc)[:600])
        finally:
            QApplication.restoreOverrideCursor()

    def _next_hotspot(self):
        """Jump to the next ambiguous stretch.

        Routed through ``_ensure_reliability`` like every other consumer of the
        analysis: this used to call ``_reliability_now``, which runs the
        perturbation ensemble *synchronously* behind a wait cursor, so the one
        command whose whole purpose is to move you around quickly was the one
        that froze the window while sixteen re-alignments ran.
        """
        if not self.aln:
            return
        self._ensure_reliability(self._show_next_hotspot)

    def _show_next_hotspot(self):
        hs = report_mod.hotspots(self.reliability.col_combined)
        if not hs:
            self.statusBar().showMessage("No ambiguous hotspots — the alignment looks clean.", 5000)
            return
        self._hotspot_idx = (getattr(self, "_hotspot_idx", -1) + 1) % len(hs)
        h = hs[self._hotspot_idx]
        self.canvas.select_nt_range(h.col_start, h.col_stop)
        self.sel_lbl.setText(
            f"hotspot {self._hotspot_idx + 1}/{len(hs)}: columns {h.col_start + 1}–{h.col_stop} "
            f"(reliability {h.mean_reliability:.2f})")

    def _compare_alignment(self, *_):
        if not self.aln:
            return
        path, _f = QFileDialog.getOpenFileName(
            self, "Compare with alignment", "",
            "Alignments (*.fasta *.fa *.fna *.faa *.aln *.phy *.sto *.nex);;All files (*)")
        if not path:
            return
        try:
            other = io.load_alignment(path)
        except Exception as exc:
            title, text = _friendly_open_error(path, exc)
            QMessageBox.critical(self, title, text)
            return
        if not (set(self.aln.ids) & set(other.ids)):
            QMessageBox.information(self, "Compare", "No shared sequence names with that "
                                   "alignment — can't compare residue placements.")
            return
        agree = disagree_mod.compare(self.aln, other)
        self.canvas.set_column_scores(agree, f"diff vs {os.path.basename(path)}")
        self.sel_lbl.setText(f"track = agreement with {os.path.basename(path)} "
                             "(red = the two alignments place residues differently)")

    # ----- manual editing (undo/redo + structural edits) -------------------- #
    def _apply_edit(self, new_aln, note: str):
        if new_aln is None:
            return
        if self.showing_truth:
            # The canvas refuses edits while the answer is up; this is the
            # backstop for any path that does not come through the canvas.
            self.sel_lbl.setText("the true alignment is read-only — Ctrl+T to go back")
            return
        self._commit(new_aln, note)               # pushes the undo entry itself

    def _undo(self):
        if not self._undo_stack:
            self.sel_lbl.setText("nothing to undo")
            return
        self._redo_stack.append(self.aln)
        self._set_alignment(self._undo_stack.pop())
        self.sel_lbl.setText("undo")

    def _redo(self):
        if not self._redo_stack:
            self.sel_lbl.setText("nothing to redo")
            return
        self._undo_stack.append(self.aln)
        self._set_alignment(self._redo_stack.pop())
        self.sel_lbl.setText("redo")

    def _cursor_cell(self):
        probe = getattr(self, "_last_probe", None)
        if not probe or not self.aln:
            return None
        seq, res = probe
        row = self.aln.rows[seq]
        cols = [c for c, ch in enumerate(row) if ch != "-"]
        return (seq, cols[res]) if res < len(cols) else None

    def _on_cursor_moved(self, row: int, col: int):
        if not self.aln:
            return
        cell = self.canvas.dm.cells[row][col]
        what = cell if cell.strip("-") else "gap"
        rows = self.canvas.selected_rows()
        who = self.aln.ids[row] if len(rows) <= 1 else f"{len(rows)} sequences"
        self.sel_lbl.setText(f"cursor · {who} · col {col + 1} ({what})")
        self._inspect_col = self.canvas.dm.nt_start[col]
        self._refresh_inspector()

    def _refresh_inspector(self):
        """Re-read the answer for the column under the cursor.

        Called when the cursor moves *and* after every edit, so that the panel
        keeps up with an alignment being worked on rather than describing the one
        it used to be.
        """
        if self.showing_truth:
            self.inspector.clear("You are looking at the true alignment itself. "
                                 "Ctrl+T returns to your alignment, where each "
                                 "column can be compared with this one.")
            return
        if not self.aln or self._inspect_col is None:
            self.inspector.clear("Click a column to see what the reference says "
                                 "about it.")
            return
        self.inspector.set_column(self.doc.explain_column(self._inspect_col))

    def _jump_to_column(self, nt_col: int):
        """Follow a link in the inspector to the column it names."""
        if not self.aln:
            return
        self._inspect_col = nt_col
        self.canvas.scroll_to_nt(nt_col)    # moves the cursor, so the panel follows
        self._refresh_inspector()

    def _on_rows_selected(self, n: int):
        if not self.aln:
            return
        if n <= 1:
            return
        self.sel_lbl.setText(f"{n} sequences selected — '-'/space slide them right, "
                             f"backspace/shift+space left")

    def _on_edit_key(self, op: str):
        """Keyboard editing acting on the *selected sequences only*. '-' / space
        slide them right (open a gap at the cursor); backspace / delete / shift+space
        slide them left (close a gap to the left). All go through undo/redo."""
        if not self.aln:
            return
        cur = self.canvas.cursor()
        if cur is None:
            self.sel_lbl.setText("click a cell first to place the edit cursor")
            return
        row, col = cur                                  # col is a *display* column
        dm = self.canvas.dm
        nt_col = dm.nt_start[col]                       # underlying nucleotide column
        unit = dm.span                                  # 1 at nt level, 3 for codon/aa
        rows = self.canvas.selected_rows()
        d = -1 if op == "slide_left" else 1
        # Moving whole codons keeps the reading frame, so keep the coding spec.
        new = editing.shift_rows(self.aln, rows, nt_col, d, unit=unit, keep_coding=unit > 1)
        if new is None:
            unitname = "codon" if unit > 1 else "gap"
            self.sel_lbl.setText(f"can't slide left — no {unitname} to the left in every "
                                 "selected sequence")
            return
        names = self.aln.ids[row] if len(rows) <= 1 else f"{len(rows)} sequences"
        self._apply_edit(new, f"slid {names} {'left' if d < 0 else 'right'} at column {col + 1}")
        self.canvas.set_cursor(row, col + d)         # cursor follows the slide
        self.canvas.set_selected_rows(rows)

    def _edit_nudge(self, direction: int):
        if not self.aln:
            return
        cell = self._cursor_cell()
        if not cell:
            self.sel_lbl.setText("click a residue first, then nudge it")
            return
        seq, col = cell
        new = editing.move_residue(self.aln, seq, col, direction)
        if new is None:
            self.sel_lbl.setText("can't nudge — no adjacent gap that way")
            return
        before = editing.column_support(self.aln, col)
        after = editing.column_support(new, col + direction)
        self._apply_edit(new, f"nudged {self.aln.ids[seq]} residue "
                              f"{'left' if direction < 0 else 'right'}")
        arrow = "↑" if after > before else ("↓" if after < before else "→")
        self.sel_lbl.setText(f"nudged · column support {before:.2f} {arrow} {after:.2f}")

    def _edit_insert_gap_col(self):
        if not self.aln:
            return
        nt = self.canvas.selection_nt()
        col = nt[0] if nt else 0
        self._apply_edit(editing.insert_gap_column(self.aln, col),
                         f"inserted gap column at {col + 1}")

    def _edit_delete_col(self):
        if not self.aln:
            return
        nt = self.canvas.selection_nt()
        col = nt[0] if nt else 0
        new = editing.delete_column(self.aln, col)
        if new is None:
            self.sel_lbl.setText(f"column {col + 1} is not all-gaps — can't delete")
            return
        self._apply_edit(new, f"deleted empty column {col + 1}")

    # ------------------------------------------------------------------ #
    # Pinned columns
    # ------------------------------------------------------------------ #
    # Pins live in ``aln.meta["anchors"]`` as a sorted list of nucleotide column
    # indices. A list rather than a set because the session writer only stores
    # JSON-safe metadata, so this shape round-trips through save/load with no
    # extra code. They are column indices into *this* alignment, so any edit that
    # changes the number of columns invalidates them; `_commit` carries them
    # across an edit and drops any that no longer exist rather than silently
    # letting them point at whatever moved into that position.

    def _pins(self) -> set:
        return set(self.aln.meta.get("anchors", ())) if self.aln else set()

    def _set_pins(self, cols) -> None:
        if not self.aln:
            return
        keep = sorted(c for c in set(cols) if 0 <= c < self.aln.length)
        if keep:
            self.aln.meta["anchors"] = keep
        else:
            self.aln.meta.pop("anchors", None)
        self.canvas.set_pins(keep)
        self._set_dirty(True)

    def _selected_nt_columns(self) -> Optional[range]:
        nt = self.canvas.selection_nt()
        return range(nt[0], nt[1]) if nt else None

    def _pin_selected(self):
        cols = self._selected_nt_columns()
        if cols is None:
            QMessageBox.information(
                self, "Pin columns",
                "Select one or more columns first — drag across the column "
                "numbers, or click a column and shift-click another.\n\n"
                "Pinned columns are held fixed by Realign around pinned columns, "
                "so you can lock the blocks you trust and let the engine re-solve "
                "only the stretches between them.")
            return
        self._set_pins(self._pins() | set(cols))
        self.sel_lbl.setText(f"{len(self._pins())} columns pinned")

    def _unpin_selected(self):
        cols = self._selected_nt_columns()
        if cols is None:
            return
        self._set_pins(self._pins() - set(cols))
        self.sel_lbl.setText(f"{len(self._pins())} columns pinned")

    def _clear_pins(self):
        self._set_pins(())
        self.sel_lbl.setText("pins cleared")

    # ------------------------------------------------------------------ #
    # Residue masks
    # ------------------------------------------------------------------ #
    # Held in ``aln.meta`` by residue (see reliability.residue_mask), so they
    # follow their residues through edits and realignment and travel in
    # sessions. Changing one is an edit: it goes on the undo stack and the
    # History, since a mask is a curation decision someone may need to retrace.

    def _selected_residues(self) -> Optional[dict]:
        """Residues under the selection: the selected sequences (or the cursor's)
        across the selected columns (or the cursor's column)."""
        rows = sorted(self.canvas.selected_rows())
        nt = self.canvas.selection_nt()
        if nt is None and self.canvas.cursor() is not None and self.canvas.dm is not None:
            c = self.canvas.cursor()[1]
            s = self.canvas.dm.nt_start[c]
            nt = (s, s + self.canvas.dm.span)
        if not rows or nt is None:
            return None
        return rel_mod.residues_in(self.aln, rows, nt[0], nt[1])

    def _mask_selected_residues(self, on: bool):
        if not self.aln:
            return
        picked = self._selected_residues()
        if picked is None:
            QMessageBox.information(
                self, "Mask residues",
                "Click a residue, or select sequences and a range of columns, "
                "first.\n\nMasked residues stay in the alignment and the column "
                "is kept; when you export, each is written as missing data (N for "
                "nucleotides, X for amino acids), so a tree program ignores it "
                "without losing the rest of the column.")
            return
        current = rel_mod.residue_mask(self.aln)
        new = rel_mod.merge_masks(current, picked) if on else rel_mod.subtract_masks(current, picked)
        changed = rel_mod.mask_size(new) - rel_mod.mask_size(current)
        if not changed:
            self.sel_lbl.setText("no change to the residue mask")
            return
        verb = "masked" if on else "unmasked"
        self._apply_edit(rel_mod.with_residue_mask(self.aln, new),
                         f"{verb} {abs(changed)} residue(s)")
        self.sel_lbl.setText(f"{verb} {abs(changed)} residue(s); "
                             f"{rel_mod.mask_size(new)} masked in all")

    def _mask_residues_below(self):
        if not self.aln:
            return
        thr = self.thr.value() / 100.0

        def apply():
            below = rel_mod.residues_below(self.aln, self.reliability.cell_combined, thr)
            current = rel_mod.residue_mask(self.aln)
            new = rel_mod.merge_masks(current, below)
            added = rel_mod.mask_size(new) - rel_mod.mask_size(current)
            if not added:
                self.sel_lbl.setText(f"no unmasked residue scores below {thr:.2f}")
                return
            self._apply_edit(rel_mod.with_residue_mask(self.aln, new),
                             f"masked {added} residue(s) with reliability < {thr:.2f}")
            self.sel_lbl.setText(f"masked {added} residue(s) below {thr:.2f}; "
                                 f"{rel_mod.mask_size(new)} masked in all")

        self._ensure_reliability(apply)

    def _clear_residue_masks(self):
        if not self.aln or not rel_mod.residue_mask(self.aln):
            return
        n = rel_mod.mask_size(rel_mod.residue_mask(self.aln))
        self._apply_edit(rel_mod.with_residue_mask(self.aln, {}), f"cleared {n} residue mask(s)")
        self.sel_lbl.setText("residue masks cleared")

    def _anchored_realign(self):
        if not self.aln:
            return
        anchors = sorted(self._pins())
        if not anchors:
            QMessageBox.information(
                self, "Realign around pinned columns",
                "No columns are pinned yet.\n\nSelect the columns you trust and "
                "pin them (Edit ▸ Edit alignment ▸ Pin selected columns, or "
                "Ctrl+P). Pinned columns are drawn with a gold bar. Realigning "
                "then holds them fixed and re-solves only the stretches between "
                "them.")
            return
        engine = self.engine_combo.currentData()
        try:
            new = editing.anchored_realign(self.aln, anchors, engine,
                                           **self._params_for(engine))
        except Exception as exc:
            QMessageBox.critical(self, "Realign failed", str(exc)[:400])
            return
        self._apply_edit(new, f"anchored realign ({len(anchors)} pinned columns)")

    def _run_agreement(self):
        """Re-align with several engines and overlay how often they reproduce
        the *current* alignment's columns — without replacing the alignment."""
        if self._busy_elsewhere("agreement"):
            self._reset_overlay("agreement")
            return
        engines = self._engines if len(self._engines) > 1 else builtin_variants()
        params_by_key = {e.key: self._params_for(e) for e in engines}
        records = self._source_records()
        alphabet = self.aln.alphabet
        base = self.aln

        def work(report, cancelled):
            return disagree_mod.compute(records, alphabet, engines, base_alignment=base,
                                        params_by_key=params_by_key, progress=report,
                                        cancelled=cancelled)

        def done(result, _secs):
            self.agreement = result.col_agreement
            self._show_track(self.agreement, "Aligner agreement")

        self._start_job("agreement", "Aligner agreement: re-aligning with each engine", work,
                        done, on_cancel=lambda: self._reset_overlay("agreement"))

    # ------------------------------------------------------------------ #
    def _on_track(self, *_):
        """The single control for the column overlay.

        Dispatch is on the track's declared kind, not on its position in the
        combo box: the list, the cheap-track rule and this function all read from
        one definition (see craic.gui.tracks), so a track can be added or
        reordered in one place.
        """
        idx = self.track_combo.currentIndex()
        if not self.aln or idx <= 0 or idx >= len(_TRACK_DEFS):
            self._show_track(None, "")
            return
        track = _TRACK_DEFS[idx]

        if track.kind == "score":
            self._show_trim_score(track.compute, track.label)
        elif track.kind == "mask":
            self._apply_binary_trim(track.compute, track.label)
        elif track.kind == "reliability":
            if track.key == "perturbation":
                self._perturbation_skipped = False     # asked for by name: offer it again
            self._ensure_reliability(lambda: self._apply_reliability_track(track))
        elif track.kind == "agreement":
            if self.agreement is not None:
                self._show_track(self.agreement, track.label)
            else:
                self._run_agreement()
        elif track.kind == "truth":
            if self.truth is None:
                self._on_track_error(
                    "No reference alignment is loaded.\n\n"
                    "Truth mode compares the current alignment with a trusted one — "
                    "a structural reference, or the known-true alignment of a "
                    "simulated dataset. Load one from the Teach menu, or generate a "
                    "dataset whose truth is known.")
                return
            self._show_track(self.truth.col_correct, "Reference correctness")
            self._report_accuracy()            # …and keep the numbers current

    def _apply_reliability_track(self, track):
        label = track.label
        if not self.reliability.n_replicates_ok:          # consistency only
            if track.key == "perturbation":
                self._reset_overlay("reliability")
                return
            if track.key == "reliability":
                label += " (consistency only)"
        self._show_track(getattr(self.reliability, track.field), label)

    def _show_track(self, scores, label):
        self.canvas.set_column_scores(scores, label)
        self._apply_mask_from_scores()
        self._update_legend()

    def _show_trim_score(self, fn, label):
        """Model-free per-column trimmer score (gap / similarity) shown as a
        mask track: the Mask slider thresholds it exactly like the other scores."""
        try:
            self._show_track(fn(self.aln), label)
        except Exception as exc:                        # pragma: no cover - defensive
            self._on_track_error(str(exc))

    def _apply_binary_trim(self, fn, label):
        """Rule-based trimmer (gappyout / strict / Gblocks) that chooses its own
        columns: show conservation for context and set the keep-mask directly."""
        try:
            keep = np.asarray(fn(self.aln), bool)
        except Exception as exc:                        # pragma: no cover - defensive
            self._on_track_error(str(exc))
            return
        self.canvas.set_column_scores(_conservation(self.aln), label)
        self.canvas.set_keep_mask(keep)
        self.sel_lbl.setText(f"{label}: keep {int(keep.sum())}/{self.aln.length} columns "
                             "(auto \u2014 move the Mask slider for threshold mode)")

    def _on_outliers(self, *_):
        """EvalMSA / OD-seq-style: flag whole sequences that agree poorly with the rest."""
        if self.aln is None or not self.outlier_chk.isChecked():
            self.canvas.set_outlier_rows(set())
            return
        try:
            _scores, flags = trim_mod.outlier_sequences(self.aln)
        except Exception as exc:                        # pragma: no cover - defensive
            self.outlier_chk.blockSignals(True)
            self.outlier_chk.setChecked(False)
            self.outlier_chk.blockSignals(False)
            QMessageBox.critical(self, "Could not flag outliers", str(exc)[-800:])
            return
        rows = {i for i, f in enumerate(flags) if f}
        self.canvas.set_outlier_rows(rows)
        self.sel_lbl.setText(f"{len(rows)} outlier sequence(s) flagged" if rows
                             else "no outlier sequences (all agree)")

    def _busy_track(self, msg):
        self.statusBar().showMessage(msg)
        self.track_combo.setEnabled(False)
        self.color_combo.setEnabled(False)

    def _idle_track(self):
        self.statusBar().clearMessage()
        self.track_combo.setEnabled(True)
        self.color_combo.setEnabled(True)

    def _on_track_error(self, msg):
        self._idle_track()
        self.track_combo.blockSignals(True)
        self.track_combo.setCurrentIndex(0)
        self.track_combo.blockSignals(False)
        self.canvas.set_column_scores(None)
        QMessageBox.critical(self, "Could not compute track", msg[-800:])

    def _on_threshold(self, *_):
        # Moving the slider IS the conservation mask: if no track is active yet,
        # turn Conservation on automatically so the slider always does something.
        if self.aln is not None and self.canvas.nt_scores is None:
            self.canvas.set_column_scores(_conservation(self.aln), "Conservation")
            self.track_combo.blockSignals(True)
            self.track_combo.setCurrentIndex(1)         # reflect it in the Track dropdown
            self.track_combo.blockSignals(False)
        self._apply_mask_from_scores()

    def _apply_mask_from_scores(self):
        """Apply the slider threshold to whatever track is currently shown."""
        thr = self.thr.value() / 100.0
        self.thr_lbl.setText(f"{thr:.2f}")
        scores = self.canvas.nt_scores
        if scores is None:
            self.canvas.set_keep_mask(None)
            return
        # One shared definition of the mask (craic.ambiguity.reliability.keep_mask),
        # so that what the user exports here is exactly what the benchmark scores.
        # Unscoreable columns are dropped, not silently kept.
        try:
            keep = rel_mod.keep_mask(scores, thr)
        except ValueError as exc:                       # no column could be scored
            self.canvas.set_keep_mask(None)
            self.sel_lbl.setText(str(exc).split(";")[0])
            return
        self.canvas.set_keep_mask(keep)
        kept = int(keep.sum())
        n_unscored = int(np.isnan(np.asarray(scores, dtype=float)).sum())
        note = f" ({n_unscored} unscored, dropped)" if n_unscored else ""
        self.sel_lbl.setText(f"mask ≥ {thr:.2f} ({self.canvas.score_label}): "
                             f"keep {kept}/{len(keep)} columns{note}")

    def _export_masked(self):
        keep = self.canvas.keep_mask
        residues = rel_mod.residue_mask(self.aln) if self.aln is not None else {}
        if self.aln is None or (keep is None and not residues):
            QMessageBox.information(self, "Nothing to mask",
                                    "Pick a Track or trimming method first "
                                    "(e.g. Reliability, Gblocks, or Conservation), "
                                    "or mask some residues.")
            return
        # Residues first, while their positions still mean something; then columns.
        masked = rel_mod.apply_residue_mask(self.aln)
        if keep is not None:
            masked = rel_mod.apply_mask(masked, np.asarray(keep, bool))
        path, _ = QFileDialog.getSaveFileName(self, "Export masked alignment",
                                              "masked.fasta", "FASTA (*.fasta)")
        if path:
            io.save_alignment(path, masked)

    # ------------------------------------------------------------------ #
    def _on_selection(self, nt0: int, nt1: int):
        self.canvas.set_probe(None)  # a new region clears a stale residue probe
        self.sandbox.set_region(nt0, nt1)
        self.posterior.set_region(nt0, nt1)
        self.sel_lbl.setText(f"selection: columns {nt0+1}–{nt1}")

    def _on_deselect(self):
        """Esc cleared the selection — forget the region in the tools too."""
        self.sandbox.clear_region()
        self.posterior.clear_region()
        self.sel_lbl.setText("selection cleared")

    # ----- named sequence groups (saved to the sidecar) --------------------- #
    def _rebuild_groups_menu(self):
        m = self.groups_menu
        m.clear()
        save = m.addAction("Save selection as group…")
        save.triggered.connect(self._save_group)
        save.setEnabled(bool(self.aln) and bool(self.canvas.selected_rows()))
        if self._groups:
            m.addSeparator()
            for name in sorted(self._groups):
                n = len(self._groups[name])
                act = m.addAction(f"Select “{name}”  ({n})")
                act.triggered.connect(lambda _c=False, nm=name: self._select_group(nm))
            m.addSeparator()
            dele = m.addAction("Delete a group…")
            dele.triggered.connect(self._delete_group)

    def _save_group(self):
        if not self.aln:
            return
        rows = sorted(self.canvas.selected_rows())
        if not rows:
            self.sel_lbl.setText("select one or more sequences first")
            return
        name, ok = QInputDialog.getText(self, "Save sequence group",
                                        "Group name (e.g. animals):")
        name = name.strip()
        if not ok or not name:
            return
        self._groups[name] = [self.aln.ids[r] for r in rows]
        self.sel_lbl.setText(f"group “{name}” = {len(rows)} sequences "
                             "(saved to the sidecar on Save)")

    def _select_group(self, name: str):
        if not self.aln or name not in self._groups:
            return
        members = set(self._groups[name])
        rows = {i for i, sid in enumerate(self.aln.ids) if sid in members}
        if not rows:
            self.sel_lbl.setText(f"group “{name}” has no sequences in this alignment")
            return
        self.canvas.set_selected_rows(rows)
        # Put the cursor on a group member (keeping the current column) so you can
        # slide straight away; clicking within the group repositions it without
        # losing the selection.
        cur = self.canvas.cursor()
        col = cur[1] if cur else 0
        first = min(rows)
        if not cur or cur[0] not in rows:
            self.canvas.set_cursor(first, col)
        missing = len(members) - len(rows)
        note = f" ({missing} not present)" if missing else ""
        self.sel_lbl.setText(f"selected group “{name}” — {len(rows)} sequences{note} · "
                             "click within the group to move the cursor")

    def _delete_group(self):
        if not self._groups:
            return
        name, ok = QInputDialog.getItem(self, "Delete group", "Group:",
                                        sorted(self._groups), 0, False)
        if ok and name in self._groups:
            del self._groups[name]
            self.sel_lbl.setText(f"deleted group “{name}”")

    # ----- column annotations (helices, active sites, …; saved to sidecar) --- #
    _ANNOT_COLORS = ["#5b8def", "#e57373", "#81c784", "#ffb74d", "#ba68c8",
                     "#4db6ac", "#f06292", "#a1887f"]

    def _rebuild_annot_menu(self):
        m = self.annot_menu
        m.clear()
        add = m.addAction("Add annotation from selection…")
        add.triggered.connect(self._add_annotation)
        add.setEnabled(bool(self.aln) and self.canvas.selection_nt() is not None)
        if self._annotations:
            m.addSeparator()
            for i, ann in enumerate(self._annotations):
                act = m.addAction(f"Select “{ann['name']}”  (cols {ann['start']+1}–{ann['stop']})")
                act.triggered.connect(lambda _c=False, k=i: self._select_annotation(k))
            m.addSeparator()
            exp = m.addAction("Export annotated columns only…")
            exp.triggered.connect(self._export_annotations)
            dele = m.addAction("Delete an annotation…")
            dele.triggered.connect(self._delete_annotation)

    def _add_annotation(self):
        if not self.aln:
            return
        rng = self.canvas.selection_nt()
        if rng is None:
            self.sel_lbl.setText("select a column range first")
            return
        name, ok = QInputDialog.getText(self, "Add annotation",
                                        "Name (e.g. helix α1, active site):")
        name = name.strip()
        if not ok or not name:
            return
        color = self._ANNOT_COLORS[len(self._annotations) % len(self._ANNOT_COLORS)]
        self._annotations.append({"name": name, "start": rng[0], "stop": rng[1],
                                  "color": color})
        self.canvas.set_annotations(self._annotations)
        self.sel_lbl.setText(f"annotation “{name}” = columns {rng[0]+1}–{rng[1]} "
                             "(saved to the sidecar on Save)")

    def _select_annotation(self, idx: int):
        if not (0 <= idx < len(self._annotations)):
            return
        ann = self._annotations[idx]
        self.canvas.select_nt_range(ann["start"], ann["stop"])   # highlights its columns
        self.sel_lbl.setText(f"annotation “{ann['name']}” — columns "
                             f"{ann['start']+1}–{ann['stop']}")

    def _delete_annotation(self):
        if not self._annotations:
            return
        labels = [f"{a['name']} (cols {a['start']+1}–{a['stop']})" for a in self._annotations]
        choice, ok = QInputDialog.getItem(self, "Delete annotation", "Annotation:",
                                          labels, 0, False)
        if ok and choice in labels:
            del self._annotations[labels.index(choice)]
            self.canvas.set_annotations(self._annotations)
            self.sel_lbl.setText("deleted annotation")

    def _export_annotations(self):
        if not self.aln or not self._annotations:
            QMessageBox.information(self, "No annotations",
                                   "Add at least one annotation first.")
            return
        keep = np.zeros(self.aln.length, dtype=bool)
        for a in self._annotations:
            keep[max(0, a["start"]):min(self.aln.length, a["stop"])] = True
        if not keep.any():
            QMessageBox.information(self, "Nothing to export", "No columns are annotated.")
            return
        masked = rel_mod.apply_mask(self.aln, keep)
        path, _ = QFileDialog.getSaveFileName(self, "Export annotated columns",
                                              "annotated.fasta", "FASTA (*.fasta)")
        if not path:
            return
        if not os.path.splitext(path)[1]:
            path += ".fasta"
        try:
            io.save_alignment(path, masked)
            self.sel_lbl.setText(f"exported {int(keep.sum())} annotated columns → {path}")
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc)[:400])

    def _busy(self, msg):
        self.statusBar().showMessage(msg)
        self.align_btn.setEnabled(False)

    def _idle(self):
        self.statusBar().clearMessage()
        self._refresh_controls()

    def _on_error(self, msg):
        self._idle()
        QMessageBox.critical(self, "Operation failed", msg[-1200:])

    def _show_about(self):
        """Who made CRAIC, where it lives, and how to cite it (with a button that
        copies the citation)."""
        box = about_box(self)
        copy_btn = box.addButton("Copy citation", QMessageBox.ActionRole)
        box.addButton(QMessageBox.Close)
        box.exec()
        if box.clickedButton() is copy_btn:
            QApplication.clipboard().setText(CITATION)
            self.statusBar().showMessage("Citation copied.", 4000)


def about_box(parent=None) -> QMessageBox:
    """The About box: version, author (linked to his website), the CRAIC website,
    the citation and the licence. Links open in the browser."""
    from .. import __author__, __version__

    box = QMessageBox(parent)
    box.setWindowTitle("About CRAIC")
    box.setTextFormat(Qt.RichText)
    box.setTextInteractionFlags(Qt.TextBrowserInteraction)
    box.setStyleSheet("QLabel{min-width: 460px;}")       # keep the links on one line
    box.setText(
        f"<h3>CRAIC {__version__}</h3>"
        f"<p>Written by <b>{__author__}</b> — "
        f"<a href='{AUTHOR_URL}'>{AUTHOR_URL}</a></p>"
        f"<p>Website and documentation: <a href='{WEBSITE}'>{WEBSITE}</a></p>"
        f"<p><b>If you use CRAIC, please cite:</b><br>{CITATION}</p>"
        f"<p>Free and open-source software under the MIT licence.</p>")
    return box


def _open_url(url: str) -> None:
    from PySide6.QtCore import QUrl
    from PySide6.QtGui import QDesktopServices

    QDesktopServices.openUrl(QUrl(url))


def _make_splash():
    """A programmatic splash / about screen (no image asset needed): program name,
    tagline, author and version, with a Continue button. Returns a frameless
    top-level QLabel whose ``continue_btn`` fires when the user clicks Continue.
    The button starts the application rather than closing it, so it says so."""
    from PySide6.QtCore import QRect, Qt
    from PySide6.QtGui import QColor, QFont, QGuiApplication, QPainter, QPixmap
    from PySide6.QtWidgets import QLabel, QPushButton

    from .. import __version__

    w, h = 540, 340
    pm = QPixmap(w, h)
    pm.fill(QColor("#14181d"))
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing, True)
    p.setRenderHint(QPainter.TextAntialiasing, True)
    p.setPen(QColor("#2a313a"))
    p.drawRect(0, 0, w - 1, h - 1)

    p.setPen(QColor("#dfe6ee"))
    f = QFont(); f.setPointSize(58); f.setBold(True); p.setFont(f)
    p.drawText(QRect(0, 46, w, 86), Qt.AlignCenter, "CRAIC")

    p.setPen(QColor("#8a97a6"))
    f2 = QFont(); f2.setPointSize(12); p.setFont(f2)
    p.drawText(QRect(20, 140, w - 40, 24), Qt.AlignCenter,
               "Conserved-Region Alignment by Iterative Convergence")

    p.setPen(QColor("#2a313a"))
    p.drawLine(170, 196, w - 170, 196)

    p.setPen(QColor("#dfe6ee"))
    f3 = QFont(); f3.setPointSize(15); p.setFont(f3)
    p.drawText(QRect(0, 208, w, 26), Qt.AlignCenter, "James McInerney")

    p.setPen(QColor("#8a97a6"))
    f4 = QFont(); f4.setPointSize(10); p.setFont(f4)
    p.drawText(QRect(0, 236, w, 20), Qt.AlignCenter,
               AUTHOR_URL.removeprefix("https://").rstrip("/"))
    p.drawText(QRect(0, 256, w, 20), Qt.AlignCenter, f"version {__version__}")
    p.end()

    splash = QLabel()
    splash.setWindowFlags(Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint)
    splash.setFixedSize(w, h)
    splash.setPixmap(pm)

    btn = QPushButton("Continue", splash)
    btn.setCursor(Qt.PointingHandCursor)
    btn.setStyleSheet(
        "QPushButton{color:#dfe6ee; background:#2a313a; border:1px solid #3a434f;"
        " border-radius:6px; padding:6px 26px;}"
        " QPushButton:hover{background:#3a434f;}")
    btn.adjustSize()
    btn.move((w - btn.width()) // 2, h - 52)
    splash.continue_btn = btn

    scr = QGuiApplication.primaryScreen().availableGeometry()
    splash.move(scr.x() + (scr.width() - w) // 2, scr.y() + (scr.height() - h) // 3)
    return splash


def _bring_to_front(win):
    """Raise and activate the window. On macOS a Qt app launched from a terminal
    often opens behind the terminal; activateWindow alone is sometimes ignored, so
    best-effort ask AppKit to bring us forward too (no-op if pyobjc is absent)."""
    import sys

    win.raise_()
    win.activateWindow()
    if sys.platform == "darwin":
        try:                                            # pragma: no cover - macOS only
            from AppKit import NSApplication
            NSApplication.sharedApplication().activateIgnoringOtherApps_(True)
        except Exception:
            pass


def _make_app_icon(side=512):
    """The app icon, drawn at runtime so there is no asset file to keep in step: a
    dark teal tile holding a grid of alignment cells over a reliability track, with
    one weak column picked out in amber, its bar dipping. Laid out on the macOS
    grid (a 10% margin round the tile) so it sits right beside other Dock icons.
    ``side`` lets the builders ask for crisp large sizes."""
    from PySide6.QtCore import QRectF, Qt
    from PySide6.QtGui import QBrush, QColor, QLinearGradient, QPainter, QPen, QPixmap

    u = side / 1024.0
    pm = QPixmap(side, side)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing, True)
    p.setPen(Qt.NoPen)

    tile = QRectF(100 * u, 92 * u, 824 * u, 824 * u)
    for i in range(14):                                  # soft drop shadow
        p.setBrush(QColor(0, 0, 0, int(10 * (1 - i / 14))))
        p.drawRoundedRect(tile.adjusted(-i * u, (14 - i) * u, i * u, (14 + i) * u),
                          (185 + i) * u, (185 + i) * u)
    grad = QLinearGradient(tile.topLeft(), tile.bottomLeft())
    grad.setColorAt(0, QColor("#1d3e48"))
    grad.setColorAt(1, QColor("#0b171c"))
    p.setBrush(QBrush(grad))
    p.drawRoundedRect(tile, 185 * u, 185 * u)
    p.setPen(QPen(QColor(255, 255, 255, 28), 3 * u))    # faint rim
    p.setBrush(Qt.NoBrush)
    p.drawRoundedRect(tile.adjusted(2 * u, 2 * u, -2 * u, -2 * u), 183 * u, 183 * u)
    p.setPen(Qt.NoPen)

    amber, teal = QColor("#f2a93b"), QColor("#45b8a6")
    light, mint = QColor("#e4edf0"), QColor("#8fd3c7")  # two quiet residue tones
    weak = 3                                             # the column the score doubts
    rows = ["ACGTTA", "ACCTTA", "AGGATA", "ACGTCA"]
    x0, cw, ch, y0 = 232 * u, 560 * u / 6, 74 * u, 262 * u
    for r, seq in enumerate(rows):
        for c, res in enumerate(seq):
            if c == weak:
                col = amber.darker(170) if r in (1, 3) else amber
            else:
                col = mint if res in "AG" else light
            p.setBrush(col)
            p.drawRoundedRect(QRectF(x0 + c * cw + 7 * u, y0 + r * (ch + 14 * u),
                                     cw - 14 * u, ch), 14 * u, 14 * u)
    base = 800 * u                                       # the reliability track
    for c, h in enumerate([0.9, 0.95, 0.85, 0.28, 0.92, 0.88]):
        p.setBrush(amber if c == weak else teal)
        p.drawRoundedRect(QRectF(x0 + c * cw + 16 * u, base - h * 150 * u,
                                 cw - 32 * u, h * 150 * u), 10 * u, 10 * u)
    p.end()
    return pm


def write_app_icon(path, side=1024):
    """Save the app icon as a PNG. The release build hands it to PyInstaller for
    the downloadable apps; needs a Qt application, and makes one if there is none
    (run with QT_QPA_PLATFORM=offscreen where there is no display)."""
    app = QApplication.instance() or QApplication([])  # noqa: F841 - Qt must exist to draw
    if not _make_app_icon(side).save(str(path), "PNG"):
        raise OSError(f"could not write the icon to {path}")


def _set_macos_dock_icon(pixmap):
    """Set the running app's Dock icon on macOS (Qt's setWindowIcon alone does not
    change the Dock tile for a non-bundled script). Needs pyobjc; no-op otherwise."""
    import sys
    if sys.platform != "darwin":
        return
    try:                                                # pragma: no cover - macOS only
        import tempfile
        from AppKit import NSApplication, NSImage
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as fh:
            path = fh.name
        pixmap.save(path, "PNG")
        img = NSImage.alloc().initWithContentsOfFile_(path)
        if img is not None:
            NSApplication.sharedApplication().setApplicationIconImage_(img)
    except Exception:
        pass


def _set_macos_app_name(name):
    """On macOS the application menu (and its Hide / Quit items) shows the process
    name — "python" / "__main__.py" — for a script that is not inside a .app bundle.
    Overriding the bundle name *before* QApplication is created makes Qt build the
    menu with our name instead. Needs pyobjc (Foundation); a no-op if it is absent."""
    import sys
    if sys.platform != "darwin":
        return
    try:                                                # pragma: no cover - macOS only
        from Foundation import NSBundle
        bundle = NSBundle.mainBundle()
        info = bundle.localizedInfoDictionary() or bundle.infoDictionary()
        if info is not None:
            info["CFBundleName"] = name
            info["CFBundleDisplayName"] = name
    except Exception:
        pass


def main():
    import sys

    from PySide6.QtCore import QTimer
    from PySide6.QtGui import QIcon
    from PySide6.QtWidgets import QApplication

    _set_macos_app_name("CRAIC")            # macOS app menu should read "CRAIC", not "python"
    app = QApplication(sys.argv)
    app.setApplicationName("CRAIC")
    app.setApplicationDisplayName("CRAIC")
    app.setStyle("Fusion")

    _icon = _make_app_icon()
    app.setWindowIcon(QIcon(_icon))
    _set_macos_dock_icon(_icon)

    try:
        splash = _make_splash()
        splash.show()
        app.processEvents()                             # paint it before we build the window
    except Exception:                                   # pragma: no cover - defensive
        splash = None

    win = CraicWindow()
    launched = {"done": False}

    def _launch():
        if launched["done"]:                            # continue + timeout both call this
            return
        launched["done"] = True
        win.show()
        _bring_to_front(win)
        if splash is not None:
            splash.close()
        # Defer the file load until the window is up and frontmost, so the
        # "not aligned" popup appears in front of it rather than behind. A file
        # named on the command line wins over a recovery offer: the user has
        # said what they want to work on.
        if len(sys.argv) > 1:
            QTimer.singleShot(0, lambda: win.open_path(sys.argv[1]))
        else:
            QTimer.singleShot(0, win.maybe_recover)

    if splash is None:
        _launch()
    else:
        splash.continue_btn.clicked.connect(_launch)    # continue early…
        QTimer.singleShot(10000, _launch)               # …or auto-continue after 10 s
    sys.exit(app.exec())
