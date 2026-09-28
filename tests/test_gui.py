"""GUI smoke test — skipped automatically if PySide6 is unavailable."""

import os
import time

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402

from craic import progressive  # noqa: E402
from craic.domain import Alignment, Alphabet, Level  # noqa: E402
from craic.gui.app import CraicWindow  # noqa: E402


def _app():
    return QApplication.instance() or QApplication([])


def test_window_constructs_and_renders():
    _app()
    win = CraicWindow()
    aln = progressive.align(
        [("a", "ATGAAAACCGCATATATTGCAAAA"), ("b", "ATGAAAACCGCATATGGGATTGCAAAA"),
         ("c", "ATGAAAACCGCGTATATTGCAAAA")],
        Alphabet.DNA,
    )
    win._set_alignment(aln)
    win.show()
    pix = win.grab()
    assert pix.width() > 0 and pix.height() > 0


def test_level_toggle_does_not_crash():
    _app()
    win = CraicWindow()
    from craic.engines import BuiltinProgressive, CodonAware

    aln = CodonAware(BuiltinProgressive()).align(
        [("a", "ATGAAAACCGCATATATTGCAAAA"), ("b", "ATGAAAACCGCATATGGGATTGCAAAA")],
        Alphabet.DNA,
    )
    win._set_alignment(aln)
    for lv in (Level.NT, Level.CODON, Level.AA):
        win.canvas.set_level(lv)
        win.canvas.grab()


def test_confidence_colour_and_probe_render():
    _app()
    from craic import progressive
    from craic.ambiguity import posterior, reliability

    win = CraicWindow()
    aln = progressive.align(
        [("a", "ACGTACGTACGT"), ("b", "ACGTAAACGTACGT"), ("c", "ACGTACGTACG")],
        Alphabet.DNA,
    )
    win._set_alignment(aln)

    rep = reliability.analyse(aln, do_perturbation=False)
    win.canvas.set_cell_scores(rep.cell_combined)
    win.canvas.set_color_mode("confidence")
    win.canvas.grab()

    win.canvas.set_probe(posterior.probe_residue(aln, 0, 2), 0, 2)
    win.canvas.grab()

    win.posterior.set_region(0, aln.length)
    for m in range(3):  # Posterior / Residual / Profile all render
        win.posterior.mode.setCurrentIndex(m)
        win.posterior.heatmap.grab()


def test_per_mode_colour():
    _app()
    from craic.gui import colors

    nt = colors.residue_color("A", Level.NT, Alphabet.DNA).name()
    cod = colors.residue_color("ATG", Level.CODON, Alphabet.DNA).name()
    aa = colors.residue_color("M", Level.AA, Alphabet.DNA).name()
    assert cod == aa and nt != cod   # codon coloured by its amino acid


def test_conservation_array():
    _app()
    from craic.gui.app import _conservation
    from craic.domain import Alignment

    aln = Alignment(["a", "b", "c"], ["AAAC", "AAAC", "AAGC"], Alphabet.DNA)
    cons = _conservation(aln)
    assert cons.shape[0] == 4
    assert cons[0] == 1.0 and abs(cons[2] - 2 / 3) < 1e-9


def test_colour_schemes_distinct():
    _app()
    from craic.gui import colors

    names = {colors.residue_color("M", Level.AA, Alphabet.PROTEIN, 1, s).name()
             for s in ("clustal", "zappo", "taylor", "hydrophobicity")}
    assert len(names) >= 3


def test_track_dropdown_selects_overlays():
    _app()
    from craic import progressive
    from craic.ambiguity import reliability

    win = CraicWindow()
    aln = progressive.align([("a", "ACGTACGT"), ("b", "ACGTAAGT"), ("c", "ACGAACGT")],
                            Alphabet.DNA)
    win._set_alignment(aln)
    win.track_combo.setCurrentText("Conservation")          # instant
    assert win.canvas._nt_scores is not None
    assert win.canvas._score_label == "Conservation"
    win.reliability = reliability.analyse(aln, do_perturbation=False)   # pre-cache
    win._perturbation_skipped = True            # as if the user chose consistency only
    win.track_combo.setCurrentText("Consistency")
    assert win.canvas._score_label == "Consistency"
    win.track_combo.setCurrentText("(none)")
    assert win.canvas._nt_scores is None


def test_new_alignment_resets_panels():
    _app()
    from craic import progressive

    win = CraicWindow()
    a = progressive.align([("x", "ACGTACGT"), ("y", "ACGTAAGT")], Alphabet.DNA)
    win._set_alignment(a)
    win._on_selection(2, 6)                       # a selection populates both panels
    assert win.sandbox.region == (2, 6)
    assert win.posterior.region == (2, 6)

    b = progressive.align([("p", "ACGTAC"), ("q", "ACGAAC")], Alphabet.DNA)
    win._set_alignment(b)                         # opening another file resets everything
    assert win.sandbox.region is None
    assert win.posterior.region is None
    assert win.sandbox.alt_list.count() == 0
    assert win.posterior.heatmap.view is None


def test_edit_invalidates_reliability_cache():
    _app()
    from craic import progressive
    from craic.ambiguity import reliability

    win = CraicWindow()
    aln = progressive.align([("a", "AAAAAA"), ("b", "AAAAAT"), ("c", "TTTTTT")], Alphabet.DNA)
    win._set_alignment(aln)
    win.reliability = reliability.analyse(aln, do_perturbation=False)   # pretend cached
    assert win.reliability is not None
    win._sort_by_similarity()                                          # an edit
    assert win.reliability is None                                     # cache discarded


def test_provenance_log_records_edits_and_splice():
    _app()
    from craic import progressive
    from craic.ambiguity import sandbox
    from craic.engines import BuiltinProgressive

    win = CraicWindow()
    aln = progressive.align(
        [("a", "ACGTACGTACGT"), ("b", "ACGTAAACGTACGT"), ("c", "ACGTACGTACG")], Alphabet.DNA)
    aln.meta["history"] = ["opened test.fasta"]
    win._set_alignment(aln)
    win._sort_by_similarity()
    assert any("sorted" in h for h in win.aln.meta["history"])
    win.sandbox.set_region(2, 8)
    win.sandbox._on_alts(sandbox.alternatives(win.aln, 2, 8, [BuiltinProgressive()]))
    win.sandbox.alt_list.setCurrentRow(0)
    win.sandbox._apply()
    assert any("realigned columns" in h for h in win.aln.meta["history"])
    assert "opened test.fasta" in win.aln.meta["history"]   # full chronological log


def test_consensus_search_and_full_image():
    _app()
    from craic import progressive

    win = CraicWindow()
    aln = progressive.align([("aa", "ACGTACGT"), ("bb", "ACGAACGT"), ("cc", "ACGTACGA")],
                            Alphabet.DNA)
    win._set_alignment(aln)
    win.canvas.set_consensus(True)
    assert win.canvas._consensus is not None and len(win.canvas._consensus) == aln.length
    assert win.canvas.search("aa") is True            # matches a sequence name
    img = win.canvas.render_full()
    assert img.width() > 0 and img.height() > 0
    win._remove_gap_columns()                          # must not crash


def test_engine_params_dialog_roundtrips_values():
    _app()
    from craic.gui.app import EngineParamsDialog

    from craic.engines import MafftEngine
    eng = MafftEngine()
    dlg = EngineParamsDialog(eng, {"op": 2.5}, None)
    vals = dlg.values()
    assert vals["op"] == 2.5                            # edited value preserved
    # the field is a typeable text box, not a spin box — typing a custom value works
    from PySide6.QtWidgets import QLineEdit
    _p, w = dlg._widgets["op"]
    assert isinstance(w, QLineEdit)
    w.setText("3.75")
    assert dlg.values()["op"] == 3.75
    dlg._restore()
    assert dlg.values()["op"] == 1.53                  # default restored


def test_a_parameter_left_to_the_tool_shows_blank_and_reads_back_as_none():
    _app()
    from PySide6.QtWidgets import QCheckBox

    from craic.engines import PrankEngine
    from craic.gui.dialogs import EngineParamsDialog, _fmt_params

    dlg = EngineParamsDialog(PrankEngine(), {})
    _p, w = dlg._widgets["gaprate"]
    assert w.text() == "" and "PRANK default" in w.placeholderText()
    vals = dlg.values()
    assert vals["gaprate"] is None and vals["gapext"] is None and vals["F"] is False
    w.setText("0.01")
    _p, box = dlg._widgets["F"]
    assert isinstance(box, QCheckBox)
    box.setChecked(True)
    vals = dlg.values()
    assert vals["gaprate"] == 0.01 and vals["F"] is True
    assert "gapext" not in _fmt_params(vals)          # History omits what PRANK chose
    dlg._restore()
    assert dlg.values()["gaprate"] is None and dlg.values()["F"] is False


def test_the_settings_dialog_shows_only_what_applies_to_the_data():
    _app()
    from craic.engines import ClustalWEngine
    from craic.gui.dialogs import EngineParamsDialog

    eng = ClustalWEngine()
    prot = EngineParamsDialog(eng, {"dnamatrix": "CLUSTALW"}, alphabet=Alphabet.PROTEIN)
    assert {"matrix", "gapdist", "nopgap", "nohgap"} <= set(prot._widgets)
    assert "dnamatrix" not in prot._widgets
    _p, w = prot._widgets["gapopen"]
    assert w.placeholderText() == "blank = ClustalW default, 10"
    assert prot.values()["dnamatrix"] == "CLUSTALW"      # hidden, but kept as left
    assert list(prot.values()) == [p.key for p in eng.parameters()]

    dna = EngineParamsDialog(eng, {}, alphabet=Alphabet.DNA)
    assert "dnamatrix" in dna._widgets
    assert not {"matrix", "gapdist", "nopgap", "nohgap"} & set(dna._widgets)
    assert dna._widgets["gapopen"][1].placeholderText() == "blank = ClustalW default, 15"

    everything = EngineParamsDialog(eng, {})              # no data loaded: show all
    assert {"matrix", "dnamatrix"} <= set(everything._widgets)
    prot._restore()
    assert prot.values()["dnamatrix"] == "IUB"             # restore reaches hidden ones


def test_the_window_asks_for_the_settings_of_the_data_the_engine_will_get(monkeypatch):
    _app()
    from craic.engines import ClustalWEngine
    from craic.gui import app as app_mod

    seen = {}

    class Recorder:
        def __init__(self, engine, values, parent=None, alphabet=None):
            seen["alphabet"] = alphabet

        def exec(self):
            return 0

    monkeypatch.setattr(app_mod, "EngineParamsDialog", Recorder)
    win = CraicWindow()
    win._set_alignment(Alignment(["a", "b"], ["ATGAAACCC", "ATGAAGCCC"], Alphabet.DNA))
    eng = ClustalWEngine()
    win.engine_combo.clear()
    win.engine_combo.addItem(eng.label, eng)
    win.codon_chk.setChecked(False)
    win._edit_engine_params()
    assert seen["alphabet"] == Alphabet.DNA
    win.codon_chk.setChecked(True)                         # align as protein
    win._edit_engine_params()
    assert seen["alphabet"] == Alphabet.PROTEIN


def test_the_settings_reach_anchored_realign_and_the_sandbox(monkeypatch):
    _app()
    from craic import editing
    from craic.engines import BuiltinProgressive

    win = CraicWindow()
    aln = progressive.align([("a", "ACGTACGTAC"), ("b", "ACGTAACGTAC"),
                             ("c", "ACGAACGTAC")], Alphabet.DNA)
    win._set_alignment(aln)
    eng = BuiltinProgressive()
    win.engine_combo.clear()
    win.engine_combo.addItem(eng.label, eng)
    win._engine_params[eng.key] = {"delta": 0.07}
    seen = {}

    def fake(aln, anchors, engine, **opts):
        seen.update(opts)
        return aln

    monkeypatch.setattr(editing, "anchored_realign", fake)
    monkeypatch.setattr(win, "_pins", lambda: [2])
    win._anchored_realign()
    assert seen["delta"] == 0.07
    assert win.sandbox.params_for(eng)["delta"] == 0.07


def test_align_shows_progress_and_recovers_on_cancel():
    _app()
    win = CraicWindow()
    win._set_alignment(progressive.align([("a", "ACGTACGT"), ("b", "ACGTAAGT"),
                                          ("c", "ACGAACGT")], Alphabet.DNA))
    win._align()
    assert win._align_progress is not None             # progress dialog appeared
    win._cancel_align()
    assert win._align_cancelled
    assert _pump_until(lambda: win._align_progress is None)   # run unwinds + closes
    assert win.align_btn.isEnabled()                   # Align is usable again


def test_params_for_overrides_defaults():
    _app()
    win = CraicWindow()

    from craic.engines import MafftEngine
    eng = MafftEngine()
    assert win._params_for(eng)["op"] == 1.53          # engine default
    win._engine_params[eng.key] = {"op": 2.0}
    assert win._params_for(eng)["op"] == 2.0           # user override wins


def _menu(win, name):
    from PySide6.QtWidgets import QToolButton

    btn = next(b for b in win.findChildren(QToolButton) if b.text() == name)
    return btn.menu()


def test_figure_export_handlers_and_hotspots(tmp_path):
    _app()
    win = CraicWindow()
    win._set_alignment(progressive.align(
        [("a", "ATGCATGCATGCAAACATGC"), ("b", "ATGCATGCATGCATGC"),
         ("c", "ATGCATGCATGCGGGCATGC")], Alphabet.DNA))
    # the report needs the reliability analysis, which runs in the background
    win._ensure_reliability(lambda: None)
    assert _pump_until(lambda: win.reliability is not None and win._job is None, 60)
    assert win.reliability.n_replicates_ok            # small data: nothing to ask about
    # the figure-export handlers (bypassing the file dialog) all produce files
    win._last_probe = (0, 4)
    base = {"seq_i": 0, "seq_j": 2, "residue": 4, "selection": False, "theme": "Light"}
    for kind, ext in [("alignment", "svg"), ("report", "png"), ("logo", "png"),
                      ("arcs", "png"), ("cooc", "png"), ("card", "png")]:
        out = str(tmp_path / f"{kind}.{ext}")
        win._render_fig({**base, "kind": kind, "fmt": ext}, out)
        assert os.path.getsize(out) > 0, kind
    # hotspot navigator selects a region (computes reliability on demand)
    win.reliability = None
    win.canvas._sel_cols = None
    win._next_hotspot()
    # either it found a hotspot (selection set) or reported none — both are valid,
    # but reliability must now be cached
    assert _pump_until(lambda: win.reliability is not None and win._job is None, 60)


def test_figure_export_dialog_picks_sequences():
    from craic.gui.app import FigureExportDialog
    _app()
    aln = progressive.align(
        [("a", "ATGCATGC"), ("b", "ATGGATGC"), ("c", "ATGCTTGC")], Alphabet.DNA)
    dlg = FigureExportDialog(aln, 1, 2, (2, 3), has_selection=True)
    kinds = [dlg.kind.itemData(i) for i in range(dlg.kind.count())]
    # paired figures expose both sequence pickers; the user's choice is honoured
    dlg.kind.setCurrentIndex(kinds.index("arcs"))
    assert dlg.seq_i.isEnabled() and dlg.seq_j.isEnabled()
    dlg.seq_i.setCurrentIndex(0); dlg.seq_j.setCurrentIndex(2)
    v = dlg.values()
    assert v["kind"] == "arcs" and v["seq_i"] == 0 and v["seq_j"] == 2
    assert v["selection"] is True and v["fmt"] in ("pdf", "svg", "png")
    # the homology card switches to a single-sequence + residue picker
    dlg.kind.setCurrentIndex(kinds.index("card"))
    assert dlg.residue.isEnabled() and not dlg.seq_j.isEnabled()


def test_manual_editing_undo_redo():
    _app()
    win = CraicWindow()
    win._set_alignment(Alignment(["a", "b"], ["AC-GT", "ACTGT"], Alphabet.DNA))
    before = list(win.aln.rows)
    win._edit_insert_gap_col()                     # widen by a gap column
    assert win.aln.length == 6 and win.aln.rows != before
    win._undo()
    assert list(win.aln.rows) == before            # back to the original
    win._redo()
    assert win.aln.length == 6                      # forward again
    # nudge the residue last clicked
    win._last_probe = (0, 2)                        # residue index 2 of seq a
    win._edit_nudge(-1)
    assert win.aln is not None                      # applied (or safely refused)


def test_keyboard_editing_acts_on_selected_rows_only():
    _app()
    win = CraicWindow()
    win._set_alignment(Alignment(["a", "b", "c"], ["ACGT-", "ACGT-", "ACGT-"], Alphabet.DNA))
    # click row b only → only that sequence slides; the others must not move
    win.canvas.set_cursor(1, 1)
    win.canvas.set_selected_rows({1})
    win._on_edit_key("slide_right")                 # '-' / space
    assert win.aln.rows[1] == "A-CGT"               # b shifted, absorbing its trailing gap
    assert win.aln.rows[0] == "ACGT-" and win.aln.rows[2] == "ACGT-"   # a, c untouched
    assert win.aln.length == 5                       # alignment did NOT widen
    win._on_edit_key("slide_left")                  # backspace / shift+space slides back
    assert win.aln.rows[1] == "ACGT-"
    # multi-select: Cmd-click adds rows; the whole set slides together
    win._set_alignment(Alignment(["a", "b", "c"], ["ACGT-", "ACGT-", "ACGTA"], Alphabet.DNA))
    win.canvas.set_cursor(0, 1)
    win.canvas.set_selected_rows({0, 1})
    win._on_edit_key("slide_right")
    assert win.aln.rows[0] == "A-CGT" and win.aln.rows[1] == "A-CGT"
    assert win.aln.rows[2] == "ACGTA"               # unselected sequence stays put
    # sliding left is refused when a selected row has no gap to its left
    win._set_alignment(Alignment(["a"], ["ACGT"], Alphabet.DNA))
    win.canvas.set_cursor(0, 2)
    win.canvas.set_selected_rows({0})
    win._on_edit_key("slide_left")
    assert win.aln.rows[0] == "ACGT"                 # unchanged, no crash


def test_editing_at_protein_level_moves_whole_codons():
    from craic.domain import CodingSpec
    _app()
    win = CraicWindow()
    win._set_alignment(Alignment(["x", "y"], ["ATGAAACCC---", "ATGAAACCCGGG"],
                                 Alphabet.DNA, coding=CodingSpec(0, 1)))
    win.canvas.set_level(Level.AA)                       # display as protein
    win.canvas.set_cursor(0, 1); win.canvas.set_selected_rows({0})
    win._on_edit_key("slide_right")                      # '-' at protein level
    assert win.aln.rows[0] == "ATG---AAACCC"             # a whole codon moved
    assert len(win.aln.rows[0]) % 3 == 0                 # frame intact
    assert win.aln.coding is not None                    # still coding
    assert win.canvas.level == Level.AA                  # view didn't snap to nt
    win._on_edit_key("slide_left")
    assert win.aln.rows[0] == "ATGAAACCC---"             # slid back


def test_shift_click_selects_row_range():
    from PySide6.QtGui import QMouseEvent
    from PySide6.QtCore import Qt, QEvent, QPointF
    _app()
    win = CraicWindow(); win.resize(1000, 400)
    win._set_alignment(Alignment([f"s{i}" for i in range(5)],
                                 ["ACGT"] * 5, Alphabet.DNA))
    cv = win.canvas
    ox, oy = cv._grid_origin(); cw = cv._col_w(); ch = cv.cell_h

    def click(r, mod=Qt.NoModifier):
        pos = QPointF(ox + 0 * cw + 3, oy + r * ch + 3)
        cv.mousePressEvent(QMouseEvent(QEvent.MouseButtonPress, pos, Qt.LeftButton,
                                       Qt.LeftButton, mod))
        cv.mouseReleaseEvent(QMouseEvent(QEvent.MouseButtonRelease, pos, Qt.LeftButton,
                                         Qt.NoButton, mod))

    click(1)                                   # anchor on s1
    click(3, Qt.ShiftModifier)                 # shift-click s3 → range s1..s3
    assert cv.selected_rows() == {1, 2, 3}
    click(0, Qt.ControlModifier)               # cmd-click adds s0
    assert cv.selected_rows() == {0, 1, 2, 3}


def test_conservation_counts_gaps_in_denominator():
    from craic.gui.app import _conservation
    _app()
    # col0 all A; col1 four T + one gap; col2 two C + three gaps; col3 all gap
    aln = Alignment(["s1", "s2", "s3", "s4", "s5"],
                    ["ATCG", "ATCG", "AT-G", "AT-G", "A--G"], Alphabet.DNA)
    cons = _conservation(aln)
    assert cons[0] == 1.0                       # fully conserved
    assert abs(cons[1] - 0.8) < 1e-9            # 80% T, 20% gap → 0.8 (not 1.0)
    assert abs(cons[2] - 0.4) < 1e-9            # gaps drag it down
    assert cons[3] == 1.0                       # all-G column


def test_mask_slider_works_without_a_track():
    _app()
    win = CraicWindow()
    win._set_alignment(Alignment(["a", "b", "c", "d"],
                                 ["ACGTACGTAC", "ACGAACTTAC", "TCGTTCGTAG", "ACGTACGTAC"],
                                 Alphabet.DNA))
    assert win.track_combo.currentIndex() == 0          # "(none)" on load
    win.thr.setValue(90)                                # move the slider only
    assert win.track_combo.currentText() == "Conservation"   # auto-engaged
    assert win.canvas._nt_keep is not None and (~win.canvas._nt_keep).any()  # some masked
    win.track_combo.setCurrentIndex(0)                  # explicitly off
    assert win.canvas._nt_keep is None                  # clears the mask


def test_column_annotations_persist_and_export(tmp_path):
    from craic import io
    _app()
    win = CraicWindow()
    win._set_alignment(Alignment(["a", "b"], ["ACGTACGTAC", "ACGTACGTAC"], Alphabet.DNA))
    # add an annotation over columns 3..6 (nt)
    win.canvas.select_nt_range(2, 6)
    win._annotations.append({"name": "helix a1", "start": 2, "stop": 6, "color": "#5b8def"})
    win.canvas.set_annotations(win._annotations)
    assert win.canvas._annot_h() > 0                      # the band is now shown
    # selecting it highlights its columns
    win.canvas.clear_selection()
    win._select_annotation(0)
    assert win.canvas.selection_nt() == (2, 6)
    # export keeps only the annotated columns
    out = str(tmp_path / "annotated.fasta")
    import craic.gui.app as app_mod
    app_mod.QFileDialog.getSaveFileName = staticmethod(lambda *a, **k: (out, ""))
    win._export_annotations()
    recs = io.read_records(out)
    assert all(len(s) == 4 for _, s in recs)             # 4 annotated columns kept
    # persists through save/open via the sidecar
    path = str(tmp_path / "aln.fasta")
    win.aln.meta["annotations"] = win._annotations
    io.save_alignment(path, win.aln)
    win.open_path(path)
    assert win._annotations and win._annotations[0]["name"] == "helix a1"


def test_click_inside_group_keeps_selection():
    from PySide6.QtGui import QMouseEvent
    from PySide6.QtCore import Qt, QEvent, QPointF
    _app()
    win = CraicWindow(); win.resize(1000, 400)
    win._set_alignment(Alignment([f"s{i}" for i in range(5)], ["ACGTACGT"] * 5, Alphabet.DNA))
    cv = win.canvas
    ox, oy = cv._grid_origin(); cw = cv._col_w(); ch = cv.cell_h

    def click(r, c):
        pos = QPointF(ox + c * cw + 3, oy + r * ch + 3)
        cv.mousePressEvent(QMouseEvent(QEvent.MouseButtonPress, pos, Qt.LeftButton,
                                       Qt.LeftButton, Qt.NoModifier))
        cv.mouseReleaseEvent(QMouseEvent(QEvent.MouseButtonRelease, pos, Qt.LeftButton,
                                         Qt.NoButton, Qt.NoModifier))

    cv.set_selected_rows({0, 1, 2})              # a group is active
    click(1, 4)                                  # plain click inside the group
    assert cv.selected_rows() == {0, 1, 2}       # group kept
    assert cv.cursor() == (1, 4)                 # cursor moved to where you clicked
    click(4, 2)                                  # plain click outside the group
    assert cv.selected_rows() == {4}             # starts a fresh single selection


def test_named_groups_save_select_persist(tmp_path):
    from craic import io
    _app()
    win = CraicWindow()
    win._set_alignment(Alignment(["cat", "dog", "fly"], ["ACGT"] * 3, Alphabet.DNA))
    # save the current selection as a group
    win.canvas.set_selected_rows({0, 1})
    win._groups["animals"] = [win.aln.ids[i] for i in (0, 1)]   # what _save_group does
    # recalling it re-selects exactly those sequences
    win.canvas.clear_selection()
    win._select_group("animals")
    assert win.canvas.selected_rows() == {0, 1}
    # groups are written to the sidecar on Save and come back on Open
    path = str(tmp_path / "a.fasta")
    win.aln.meta["groups"] = win._groups
    io.save_alignment(path, win.aln)
    win.open_path(path)
    assert win._groups == {"animals": ["cat", "dog"]}


def test_escape_clears_selection():
    _app()
    win = CraicWindow()
    win._set_alignment(Alignment(["a", "b"], ["ACGTACGT", "ACGTACGT"], Alphabet.DNA))
    win.canvas.select_all()
    win.canvas.set_selected_rows({0})
    assert win.canvas._sel_cols is not None and win.canvas.selected_rows() == {0}
    win.sandbox.set_region(0, 8)
    win.canvas.clear_selection()                 # the Esc key
    assert win.canvas._sel_cols is None
    assert win.canvas._sel_rows == set()
    assert win.sandbox.region is None            # the tools forgot the region too
    assert win.posterior.region is None


def test_frame_break_warning_fires_once(monkeypatch):
    from craic.domain import CodingSpec
    from craic.gui import app as app_mod
    _app()
    win = CraicWindow()
    # gaps not on codon boundaries → partial-gap codons that translate to X
    win._set_alignment(Alignment(["a", "b"], ["ATG-AAACCC", "ATGAA-ACCC"],
                                 Alphabet.DNA, coding=CodingSpec(0, 1)))
    calls = []
    monkeypatch.setattr(app_mod.QMessageBox, "warning",
                        lambda *a, **k: calls.append(a))
    win.canvas.set_level(Level.AA)
    win._maybe_warn_frame_breaks()
    assert len(calls) == 1 and win._warned_frame_breaks
    win._maybe_warn_frame_breaks()               # doesn't nag again
    assert len(calls) == 1


def test_font_zoom_scales_cells_and_text():
    _app()
    win = CraicWindow()
    win._set_alignment(Alignment(["a", "b"], ["ACGT", "ACGT"], Alphabet.DNA))
    cv = win.canvas
    h0, f0 = cv.cell_h, cv._letter.pixelSize()
    cv.zoom(3)                                   # Ctrl/Cmd + '+' (or the View menu)
    assert cv.cell_h > h0 and cv._letter.pixelSize() > f0   # cells AND text grow
    cv.zoom(-3)
    assert cv.cell_h == h0
    cv.set_cell_size(40, 44)                      # clamps at the top end
    assert cv.cell_h == 44
    cv.reset_zoom()
    assert (cv.cell_w, cv.cell_h) == (15, 17)


def test_mouse_drag_slides_selected_rows():
    from PySide6.QtGui import QMouseEvent
    from PySide6.QtCore import Qt, QEvent, QPointF
    _app()
    win = CraicWindow(); win.resize(1000, 400)
    win._set_alignment(Alignment(["a", "b", "c"],
                                 ["ACGTACGT", "--ACGTAC", "ACGTACGT"], Alphabet.DNA))
    cv = win.canvas
    cv.viewport().resize(800, 340)
    ox, oy = cv._grid_origin(); cw = cv._col_w(); ch = cv.cell_h

    def ev(kind, r, c, btn=Qt.LeftButton, btns=Qt.LeftButton):
        pos = QPointF(ox + c * cw + 3, oy + r * ch + 3)
        return QMouseEvent(kind, pos, btn, btns, Qt.NoModifier)

    # grab the (already-selected) sequence b and drag it right two columns
    cv.set_selected_rows({1}); cv.set_cursor(1, 4)
    cv.mousePressEvent(ev(QEvent.MouseButtonPress, 1, 4))
    cv.mouseMoveEvent(ev(QEvent.MouseMove, 1, 6, btn=Qt.NoButton))
    cv.mouseReleaseEvent(ev(QEvent.MouseButtonRelease, 1, 6, btn=Qt.NoButton))
    # b's residues slid right (gaps moved), and no residue of a or c moved — they
    # only gain harmless trailing padding to keep the block rectangular
    assert win.aln.rows[1].replace("-", "") == "ACGTAC"
    assert win.aln.rows[1] != "--ACGTAC"
    assert win.aln.rows[0].replace("-", "") == "ACGTACGT"
    assert win.aln.rows[2].replace("-", "") == "ACGTACGT"
    assert win.aln.rows[0].rstrip("-") == "ACGTACGT"             # body unchanged
    # drag it back left to where it started
    cv.mousePressEvent(ev(QEvent.MouseButtonPress, 1, 6))
    cv.mouseMoveEvent(ev(QEvent.MouseMove, 1, 4, btn=Qt.NoButton))
    cv.mouseReleaseEvent(ev(QEvent.MouseButtonRelease, 1, 4, btn=Qt.NoButton))
    assert win.aln.rows[1].replace("-", "") == "ACGTAC"          # still intact


def test_view_menu_toggles_ambiguity_dock():
    from PySide6.QtWidgets import QToolButton

    _app()
    win = CraicWindow(); win.show()
    assert not win.dock.isVisible()                    # hidden by default
    btn = next(b for b in win.findChildren(QToolButton) if b.text() == "View")
    toggle = next(a for a in btn.menu().actions() if "Ambiguity" in a.text())
    toggle.trigger()
    assert win.dock.isVisible()                        # View menu brings it up
    toggle.trigger()
    assert not win.dock.isVisible()                    # and hides it again


def test_file_and_edit_menus_group_actions():
    _app()
    win = CraicWindow()
    file_labels = [a.text() for a in _menu(win, "File").actions() if not a.isSeparator()]
    edit_labels = [a.text() for a in _menu(win, "Edit").actions() if not a.isSeparator()]
    assert file_labels == ["Open…", "Open session…",
                           "Save session", "Save session as…", "Save alignment as…",
                           "Export image…", "Export masked…",
                           "Export figure…", "Compare with alignment…", "History…"]
    assert edit_labels == ["Copy", "Select all", "Sort by similarity", "Remove gap columns",
                           "Go to next ambiguous region", "Edit alignment", "Sequence groups",
                           "Annotations"]
    sc = {a.text(): a.shortcut().toString() for a in _menu(win, "File").actions()}
    # Ctrl+S saves the *session*: it is the one that loses nothing, so it is what
    # a reflex should reach. Exporting the alignment alone is Ctrl+Shift+S.
    assert sc["Open…"] == "Ctrl+O"
    assert sc["Save session"] == "Ctrl+S"
    assert sc["Save alignment as…"] == "Ctrl+Shift+S"


def test_document_actions_disabled_without_alignment():
    _app()
    win = CraicWindow()
    acts = {a.text(): a for a in _menu(win, "File").actions() + _menu(win, "Edit").actions()}
    assert acts["Save alignment as…"].isEnabled() is False     # nothing loaded yet
    assert acts["Save session"].isEnabled() is False
    win._set_alignment(progressive.align([("a", "ACGT"), ("b", "ACGT")], Alphabet.DNA))
    assert acts["Save alignment as…"].isEnabled() is True
    assert acts["Save session"].isEnabled() is True
    assert acts["Sort by similarity"].isEnabled() is True


def test_copy_submenu_offers_formats_and_copies():
    _app()
    from PySide6.QtWidgets import QApplication

    win = CraicWindow()
    win._set_alignment(progressive.align([("a", "ACGTACGT"), ("b", "ACGAACGT")], Alphabet.DNA))
    copy_menu = next(a.menu() for a in _menu(win, "Edit").actions() if a.menu() is not None)
    labels = [a.text() for a in copy_menu.actions()]
    assert {"FASTA", "Clustal", "MEGA"} <= set(labels)
    sc = {a.text(): a.shortcut().toString() for a in copy_menu.actions()}
    assert sc["FASTA"] == "Ctrl+C"                     # quick-copy still FASTA
    win._copy_as("mega")
    assert "#MEGA" in QApplication.clipboard().text()
    win._copy_as("fasta")
    assert QApplication.clipboard().text().startswith(">a")


def test_select_all_covers_whole_alignment():
    _app()
    win = CraicWindow()
    win._set_alignment(progressive.align([("a", "ACGTACGTAC"), ("b", "ACGAACGTAC")], Alphabet.DNA))
    cv = win.canvas
    cv.select_all()
    assert cv._sel_cols == (0, cv._dm.n_cols)
    sa = next(a for a in _menu(win, "Edit").actions() if a.text() == "Select all")
    assert sa.shortcut().toString() == "Ctrl+A"


def test_drag_past_edge_scrolls_and_extends_selection():
    _app()
    from PySide6.QtWidgets import QApplication

    win = CraicWindow(); win.resize(1200, 400)          # wide enough that the grid is visible
    row = "ACGT" * 80                                   # 320 columns → needs scrolling
    win._set_alignment(Alignment(["a", "b"], [row, row], Alphabet.DNA))
    win.show(); QApplication.processEvents()
    cv = win.canvas
    assert cv.viewport().width() > cv.name_w            # sanity: there is a visible grid
    cv._anchor = 0; cv._moved = False
    cv._drag_x = cv.viewport().width() - 1              # cursor pinned to the right edge
    hbar = cv.horizontalScrollBar(); start = hbar.value()
    cv._update_autoscroll(cv._drag_x)
    assert cv._autoscroll_dir == 1 and cv._autoscroll.isActive()
    for _ in range(5):
        cv._autoscroll_tick()
    cv._autoscroll.stop()
    assert hbar.value() > start                         # the view scrolled
    assert cv._sel_cols and cv._sel_cols[1] > 1         # and the selection grew


def test_posterior_colormaps_distinct_and_safe():
    from craic.gui import colors

    assert {"Viridis", "Greyscale"} <= set(colors.COLORMAP_NAMES)
    a = colors.colormap_color(0.5, "Viridis").name()
    b = colors.colormap_color(0.5, "Magma").name()
    grey = colors.colormap_color(0.6, "Greyscale")
    assert a != b                                       # different maps → different colours
    assert grey.red() == grey.green() == grey.blue()    # greyscale really is grey
    # Magma (pink/salmon) and Inferno (orange/gold) must be visibly different in
    # the mid/high range, not just at the shared endpoints.
    for v in (0.625, 0.75, 0.875):
        m = colors.colormap_color(v, "Magma")
        f = colors.colormap_color(v, "Inferno")
        assert abs(m.green() - f.green()) + abs(m.blue() - f.blue()) > 40
    for v in (0.0, 1.0, float("nan")):                  # endpoints + nan are safe
        assert colors.colormap_color(v, "Turbo").isValid()


def test_posterior_panel_colormap_switch():
    _app()
    from craic.gui import colors
    from craic.gui.panels import PosteriorPanel

    pp = PosteriorPanel()
    assert [pp.cmap.itemText(i) for i in range(pp.cmap.count())] == colors.COLORMAP_NAMES
    assert pp.heatmap.cmap == "Viridis"                 # colourful by default
    pp.cmap.setCurrentText("Magma")
    assert pp.heatmap.cmap == "Magma"                   # dropdown drives the heatmap


def test_posterior_panel_caps_huge_region():
    _app()
    from craic.gui.panels import MAX_POSTERIOR_COLS, PosteriorPanel

    row = "ACGT" * 70                                   # 280 columns
    pp = PosteriorPanel()
    pp.set_alignment(Alignment(["a", "b", "c"], [row, row, row], Alphabet.DNA))
    pp.set_region(0, 280)                               # "Select all" — whole alignment
    assert pp.heatmap.view is None                      # refused: no giant matrix built/painted
    assert "too large" in pp.info.text().lower()
    pp.set_region(0, MAX_POSTERIOR_COLS // 2)           # a local region still works
    assert pp.heatmap.view is not None


def test_plain_click_probes_residue():
    _app()
    win = CraicWindow()
    win._set_alignment(progressive.align([("a", "ACGTACGT"), ("b", "ACGAACGT")], Alphabet.DNA))
    cv = win.canvas
    probed = []
    cv.residueClicked.connect(lambda *a: probed.append(a))
    # a plain click (no drag) probes the clicked residue
    cv._anchor, cv._press, cv._moved = 2, (2, 0), False
    cv.mouseReleaseEvent(None)
    assert probed
    # a drag emits a selection, not a probe
    probed.clear()
    sel = []
    cv.selectionChanged.connect(lambda a, b: sel.append((a, b)))
    cv._anchor, cv._press, cv._moved = 1, (1, 0), True
    cv._sel_cols = (1, 5)
    cv.mouseReleaseEvent(None)
    assert sel and not probed


def test_figures_export_every_kind(tmp_path):
    _app()
    from craic import figures
    from craic.ambiguity import ensemble, reliability as rel_mod
    from craic.engines import builtin_variants

    aln = progressive.align([("a", "ACGTACGTACGT"), ("b", "ACGTAAACGTAC"),
                             ("c", "ACGTACGTACGT")], Alphabet.DNA)
    rep = rel_mod.analyse(aln, do_perturbation=False)
    opt = figures.FigureOptions(reliability=rep.col_combined, title="t")

    checks = [
        ("align.svg", lambda p: figures.export_alignment(aln, p, opt=opt)),
        ("align.pdf", lambda p: figures.export_alignment(aln, p, opt=opt)),
        ("card.png", lambda p: figures.homology_card(aln, 0, 2, p)),
        ("logo.png", lambda p: figures.uncertainty_logo(aln, 0, aln.length, p,
                                                        reliability=rep.col_combined)),
        ("arcs.png", lambda p: figures.homology_arcs(aln, 0, 1, 0, aln.length, p)),
        ("report.png", lambda p: figures.confidence_report(aln, rep.col_combined, p)),
        ("view.html", lambda p: figures.export_html(aln, p, reliability=rep.col_combined)),
    ]
    for name, fn in checks:
        out = str(tmp_path / name)
        fn(out)
        assert os.path.getsize(out) > 0, name

    co, ci, cj, _n = ensemble.cooccurrence(aln, 0, 1, 0, aln.length, builtin_variants())
    out = str(tmp_path / "cooc.svg")
    figures.matrix_heatmap(co, list(ci), list(cj), out)
    assert os.path.getsize(out) > 0


def _pump_until(predicate, seconds=5):
    from PySide6.QtWidgets import QApplication

    deadline = time.time() + seconds
    while not predicate() and time.time() < deadline:
        QApplication.processEvents()
        time.sleep(0.01)
    return predicate()


def test_run_cancellable_reports_progress_and_completes():
    _app()
    from craic.gui.workers import run_cancellable

    state = {"done": None, "progress": []}

    def work(report, cancelled):
        report(1, 2, "step 1")
        report(2, 2, "step 2")
        return "ok"

    run_cancellable(work, lambda r: state.update(done=r),
                    on_progress=lambda d, t, lbl: state["progress"].append((d, t, lbl)))
    assert _pump_until(lambda: state["done"] == "ok")
    assert (2, 2, "step 2") in state["progress"]


def test_sandbox_generate_shows_progress_then_cancels():
    _app()
    win = CraicWindow()
    aln = progressive.align([("a", "ACGTACGTACGT"), ("b", "ACGTAAACGTACGT"),
                             ("c", "ACGTACGTACG")], Alphabet.DNA)
    win._set_alignment(aln)
    sb = win.sandbox
    sb.set_region(0, aln.length)
    sb._generate()
    assert sb._progress is not None          # the progress dialog appeared
    sb._cancel_generate()
    assert sb._cancelled                      # cancel requested
    assert _pump_until(lambda: sb._progress is None)   # run unwinds and closes


def test_score_trimmers_are_mask_tracks():
    """Gap-fraction and similarity trimmers behave like any score track:
    they populate the band and the slider turns them into a keep-mask."""
    _app()
    from craic import progressive

    win = CraicWindow()
    aln = progressive.align([("a", "ACGTACGTAC"), ("b", "ACGT--GTAC"), ("c", "ACGTACGTAG")],
                            Alphabet.DNA)
    win._set_alignment(aln)
    for name in ("Gap fraction (MSA_trimmer)", "Similarity (trimAl)"):
        win.track_combo.setCurrentText(name)
        assert win.canvas._score_label == name
        assert win.canvas.nt_scores is not None
        assert win.canvas.keep_mask is not None
        assert len(win.canvas.keep_mask) == aln.length


def test_rule_based_trimmers_set_keep_mask():
    """gappyout / strict / Gblocks choose their own columns -> a bool keep-mask."""
    _app()
    from craic import progressive

    win = CraicWindow()
    aln = progressive.align([("a", "ACGTACGTAC"), ("b", "ACGT--GTAC"), ("c", "ACGTACGTAG")],
                            Alphabet.DNA)
    win._set_alignment(aln)
    for name in ("trimAl gappyout", "trimAl strict", "Gblocks"):
        win.track_combo.setCurrentText(name)
        km = win.canvas.keep_mask
        assert km is not None and len(km) == aln.length
        assert km.dtype == bool


def test_outlier_checkbox_flags_and_clears_rows():
    """The outlier toggle runs EvalMSA-style detection and updates the canvas;
    unchecking clears the highlight."""
    _app()
    from craic import progressive

    win = CraicWindow()
    aln = progressive.align(
        [("a", "ACGTACGTACGTACGT"), ("b", "ACGTACGTACGTACGT"),
         ("c", "ACGTACGTACGTACGT"), ("odd", "TTTTGGGGCCCCAAAA")],
        Alphabet.DNA)
    win._set_alignment(aln)
    win.outlier_chk.setChecked(True)
    assert isinstance(win.canvas._outlier_rows, set)
    win.outlier_chk.setChecked(False)
    assert win.canvas._outlier_rows == set()


def test_friendly_open_error_is_concise():
    """A parser error that embeds a whole sequence must become a short popup."""
    from craic.gui.app import _friendly_open_error

    huge = "0100" * 2000
    exc = ValueError(f"Taxon AEDAE: Illegal character 1 in sequence {huge}")
    title, text = _friendly_open_error("/data/gene_content.nex", exc)
    assert "gene_content.nex" in title
    assert huge not in text                    # the giant sequence is gone
    assert len(text) < 400
    assert "AEDAE" in text and "presence-absence" in text


def test_an_edit_is_one_undo_step():
    """Hand edits used to push the previous alignment twice, so the second Undo
    after one edit appeared to do nothing."""
    _app()
    from craic import editing

    win = CraicWindow()
    win._set_alignment(Alignment(["a", "b"], ["AC-GT", "ACTGT"], Alphabet.DNA))
    before = win.aln
    win._apply_edit(editing.insert_gap_column(win.aln, 1), "gap")
    assert len(win._undo_stack) == 1
    win._undo()
    assert win.aln is before and not win._undo_stack


def test_masking_residues_is_an_undoable_edit_that_survives_realigning(tmp_path, monkeypatch):
    _app()
    from PySide6.QtWidgets import QFileDialog

    from craic import io
    from craic.ambiguity import reliability as rel_mod
    from craic.engines import BuiltinProgressive

    win = CraicWindow()
    win._set_alignment(progressive.align(
        [("a", "ATGAAAACCGCATATATTGCAAAA"), ("b", "ATGAAAACCGCATATGGGATTGCAAAA"),
         ("c", "ATGAAAACCGCGTATATTGCAAAA")], Alphabet.DNA))
    win.canvas.set_selected_rows({1})
    win.canvas.select_nt_range(12, 18)
    win._mask_selected_residues(True)
    mask = rel_mod.residue_mask(win.aln)
    assert list(mask) == ["b"] and len(mask["b"]) == 6
    assert win.aln.meta["history"][-1] == "masked 6 residue(s)"
    assert win.canvas._cell_masked(1, 13)

    records = [(i, r.replace("-", "")) for i, r in zip(win.aln.ids, win.aln.rows)]
    win._commit(BuiltinProgressive().align(records, Alphabet.DNA, effort="min"), "aligned")
    assert rel_mod.residue_mask(win.aln) == mask          # the residues, not the columns

    out = tmp_path / "masked.fasta"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(out), ""))
    win._export_masked()
    back = io.load_alignment(str(out))
    assert back.length == win.aln.length                  # no column dropped
    assert back.rows[1].replace("-", "")[12:18] == "NNNNNN"

    win._undo()                                           # back past the realignment
    win._undo()                                           # and past the mask
    assert not rel_mod.residue_mask(win.aln)


def test_help_menu_credits_author_and_citation():
    from craic import AUTHOR_URL, CITATION, WEBSITE, __author__
    from craic.gui.app import about_box

    _app()
    win = CraicWindow()
    texts = [a.text() for a in _menu(win, "Help").actions()]
    assert texts == ["CRAIC website", "About CRAIC and how to cite it…"]
    box = about_box(win)
    for s in (__author__, AUTHOR_URL, WEBSITE, CITATION):
        assert s in box.text()


def test_write_app_icon(tmp_path):
    from PySide6.QtGui import QImage

    from craic.gui.app import write_app_icon

    _app()
    out = tmp_path / "icon.png"
    write_app_icon(out, 256)
    img = QImage(str(out))
    assert (img.width(), img.height()) == (256, 256)
    assert img.pixelColor(0, 0).alpha() == 0             # transparent round the tile


# --------------------------------------------------------------------------- #
# Long analyses: progress, the up-front question, cancel (0.5.11)
# --------------------------------------------------------------------------- #

def _three_seq_window():
    win = CraicWindow()
    win._set_alignment(progressive.align(
        [("a", "ATGCATGCATGCAAACATGC"), ("b", "ATGCATGCATGCATGC"),
         ("c", "ATGCATGCATGCGGGCATGC")], Alphabet.DNA))
    return win


def test_a_long_perturbation_run_is_offered_not_imposed():
    _app()
    win = _three_seq_window()
    win._ASK_ABOVE_S = 0                        # every run counts as long here
    asked = []
    win._confirm_perturbation = lambda secs: asked.append(secs) or False
    win.track_combo.setCurrentText("Reliability")
    assert _pump_until(lambda: win._job is None and win.reliability is not None, 30)
    assert asked and win.reliability.n_replicates_ok == 0      # consistency only
    assert win.canvas._score_label == "Reliability (consistency only)"

    # asking for the perturbation track by name offers the run again
    win._confirm_perturbation = lambda secs: True
    win.track_combo.setCurrentText("Perturbation")
    assert _pump_until(lambda: win._job is None and win.reliability.n_replicates_ok > 0, 60)
    assert win.canvas._score_label == "Perturbation"


def test_cancelling_the_analysis_leaves_the_window_as_it_was():
    _app()
    win = _three_seq_window()
    win.track_combo.setCurrentText("Reliability")
    assert win._job is not None
    win._cancel_job(win._job)
    assert _pump_until(lambda: win._job is None, 30)
    assert win.reliability is None and not win._reliability_waiters
    assert win.track_combo.currentIndex() == 0 and win.track_combo.isEnabled()


def test_an_edit_drops_the_running_analysis():
    _app()
    win = _three_seq_window()
    win.track_combo.setCurrentText("Reliability")
    assert win._job is not None
    win._set_alignment(progressive.align([("x", "ACGTACGT"), ("y", "ACGTAAGT")], Alphabet.DNA))
    assert win._job is None and win.reliability is None
    _pump_until(lambda: False, 1)               # the old job's result arrives and is dropped
    assert win.reliability is None


def test_progress_dialog_says_how_long_is_left():
    _app()
    from craic.gui.dialogs import JobProgress, duration

    assert duration(20) == "less than a minute"
    assert duration(40 * 60) == "40 min"
    assert duration(3 * 3600) == "3.0 h"
    dlg = JobProgress("CRAIC", "Reliability")
    dlg._t0 -= 60                                # a minute in, a quarter done
    dlg.update_progress(25, 100)
    assert "About 3 min left" in dlg.labelText()
    assert dlg.value() == 250
