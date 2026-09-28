"""Tests for the v0.3 additions: the headless CLI, truth mode, the dataset
generator, and the sparse consistency transformation."""

import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from craic import evaluate, progressive, simulate
from craic.domain import Alignment, Alphabet, IdMismatch

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _sim(**kw):
    kw.setdefault("taxa", 5)
    kw.setdefault("root_len", 120)
    kw.setdefault("seed", 1)
    d = simulate.simulate(**kw)
    truth = Alignment(list(d["names"]), list(d["true_rows"]), Alphabet.DNA)
    return d, truth


# --------------------------------------------------------------------------- #
# evaluate: accuracy against a reference
# --------------------------------------------------------------------------- #

def test_an_alignment_scores_perfectly_against_itself():
    _d, truth = _sim()
    acc = evaluate.compare_to_reference(truth, truth)
    assert acc.sp == 1.0 and acc.tc == 1.0 and acc.precision == 1.0
    assert np.nanmin(acc.col_correct) == 1.0


def test_reference_rows_are_matched_by_id_not_by_position():
    """A reference that lists the same sequences in a different order is the
    normal case, not an error."""
    _d, truth = _sim()
    order = list(range(truth.n_seqs))[::-1]
    shuffled = Alignment([truth.ids[i] for i in order],
                         [truth.rows[i] for i in order], Alphabet.DNA)
    assert evaluate.compare_to_reference(truth, shuffled).sp == 1.0


def test_a_reference_for_different_sequences_is_an_error_not_a_low_score():
    """Otherwise a reference for the wrong data reads as an alignment that is
    entirely wrong, which is a far more damaging way to fail."""
    _d, truth = _sim(seed=1)
    _d2, other = _sim(seed=2)
    other = Alignment(list(truth.ids), list(other.rows[:truth.n_seqs]), Alphabet.DNA)
    with pytest.raises(IdMismatch):
        evaluate.compare_to_reference(truth, other)


def test_misaligned_columns_are_marked_and_reliability_auc_is_defined():
    d, truth = _sim(taxa=6, root_len=150, seed=7, indel_rate=3.0, bmax=0.4)
    from craic.engines import BuiltinProgressive
    from craic.ambiguity import reliability as rel

    aln = BuiltinProgressive().align(d["seqs"], Alphabet.DNA)
    acc = evaluate.compare_to_reference(aln, truth)
    assert 0.0 < acc.sp < 1.0                       # a genuinely hard dataset
    assert acc.col_correct.size == aln.length
    rep = rel.analyse(aln, do_perturbation=False)
    a = acc.reliability_auc(rep.col_consistency)
    assert a == a and 0.0 <= a <= 1.0


# --------------------------------------------------------------------------- #
# The command line
# --------------------------------------------------------------------------- #

def _craic(*args, cwd):
    env = dict(os.environ, PYTHONPATH=REPO)
    return subprocess.run([sys.executable, "-m", "craic", *args],
                          cwd=cwd, env=env, capture_output=True, text=True, timeout=600)


def test_cli_round_trip_simulate_align_score_mask(tmp_path):
    p = _craic("simulate", "--taxa", "5", "--length", "100", "--seed", "2",
               "-o", "truth.fasta", "--unaligned", "seqs.fasta", cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    assert (tmp_path / "truth.fasta").exists() and (tmp_path / "seqs.fasta").exists()

    p = _craic("align", "seqs.fasta", "-o", "aln.fasta", cwd=tmp_path)
    assert p.returncode == 0, p.stderr

    p = _craic("score", "aln.fasta", "--fast", "-o", "scores.tsv",
               "--reference", "truth.fasta", cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    assert "SP=" in p.stderr and "reliability_AUC=" in p.stderr
    head, *rows = (tmp_path / "scores.tsv").read_text().splitlines()
    assert head.split("\t")[0] == "column"
    assert "correct_fraction" in head          # the reference adds truth to the table
    assert len(rows) > 0

    p = _craic("mask", "aln.fasta", "--fast", "--threshold", "0.5",
               "-o", "masked.fasta", cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    assert "kept" in p.stderr


def test_cli_writes_to_stdout_when_no_output_file(tmp_path):
    _craic("simulate", "--taxa", "4", "--length", "80", "-o", "truth.fasta", cwd=tmp_path)
    p = _craic("trim", "truth.fasta", "--method", "gappyout", cwd=tmp_path)
    assert p.returncode == 0, p.stderr
    assert p.stdout.startswith(">")            # FASTA on stdout, progress on stderr


def test_cli_refuses_a_reference_for_the_wrong_sequences(tmp_path):
    _craic("simulate", "--taxa", "4", "--length", "80", "--seed", "1",
           "-o", "a.fasta", cwd=tmp_path)
    _craic("simulate", "--taxa", "4", "--length", "80", "--seed", "9",
           "-o", "b.fasta", cwd=tmp_path)
    p = _craic("score", "a.fasta", "--fast", "--reference", "b.fasta",
               "-o", "-", cwd=tmp_path)
    assert p.returncode != 0
    assert "craic:" in p.stderr


def test_cli_reports_the_backend(tmp_path):
    p = _craic("--core", cwd=tmp_path)
    assert p.returncode == 0
    assert p.stdout.strip() in ("rust", "numpy")


def test_headless_subcommands_never_import_qt(tmp_path):
    """The reason the CLI exists: it has to run where there is no Qt at all."""
    _craic("simulate", "--taxa", "4", "--length", "80", "-o", "t.fasta", cwd=tmp_path)
    script = textwrap.dedent("""
        import sys
        class Block:
            def find_module(self, name, path=None):
                if name == "PySide6" or name.startswith("PySide6."):
                    return self
            def load_module(self, name):
                raise ImportError("no Qt here")
        sys.meta_path.insert(0, Block())
        from craic.cli import main
        try:
            main(["score", "t.fasta", "--fast", "-o", "out.tsv"])
        except SystemExit:
            pass
        leaked = [m for m in sys.modules if m.startswith("PySide6")]
        assert not leaked, leaked
        print("clean")
    """)
    env = dict(os.environ, PYTHONPATH=REPO)
    p = subprocess.run([sys.executable, "-c", script], cwd=tmp_path, env=env,
                       capture_output=True, text=True, timeout=600)
    assert p.returncode == 0, p.stderr
    assert "clean" in p.stdout


def test_bare_invocation_is_still_the_workbench():
    """`craic` and `craic FILE` must keep launching the GUI: the subcommands are
    additive, not a replacement."""
    from craic import cli

    launched = []
    original = cli._launch_gui
    cli._launch_gui = lambda argv: launched.append(list(argv))
    try:
        cli.main([])
        cli.main(["some.fasta"])
    finally:
        cli._launch_gui = original
    assert launched == [[], ["some.fasta"]]


# --------------------------------------------------------------------------- #
# The sparse consistency transformation
# --------------------------------------------------------------------------- #

@pytest.mark.skipif(not progressive.sparse_available(), reason="SciPy not installed")
def test_sparse_and_dense_transforms_agree_closely():
    d, _truth = _sim(taxa=5, root_len=150, seed=3)
    seqs = [s for _, s in d["seqs"]]
    model, delta, eps = progressive.estimate_params(d["true_rows"], Alphabet.DNA)
    raw = progressive._all_pairs_posteriors(seqs, model, delta, eps)
    n = len(seqs)
    dense = progressive.consistency_transform(dict(raw), n, 2, sparse=False)
    sparse = progressive.consistency_transform(dict(raw), n, 2, sparse=True)
    worst = max(np.abs(progressive._post(dense, x, y) - progressive._post(sparse, x, y)).max()
                for x in range(n) for y in range(x + 1, n))
    # Thresholding at 0.01 is a deliberate approximation; it must stay small.
    assert worst < 0.1


@pytest.mark.skipif(not progressive.sparse_available(), reason="SciPy not installed")
def test_sparse_transform_stores_far_less():
    d, _truth = _sim(taxa=5, root_len=200, seed=3)
    seqs = [s for _, s in d["seqs"]]
    model, delta, eps = progressive.estimate_params(d["true_rows"], Alphabet.DNA)
    raw = progressive._all_pairs_posteriors(seqs, model, delta, eps)
    n = len(seqs)
    dense = progressive.consistency_transform(dict(raw), n, 1, sparse=False)
    sparse = progressive.consistency_transform(dict(raw), n, 1, sparse=True)
    db = sum(v.nbytes for v in dense.values())
    sb = sum(v.data.nbytes + v.indices.nbytes + v.indptr.nbytes for v in sparse.values())
    assert sb < 0.5 * db            # in practice a few per cent


@pytest.mark.skipif(not progressive.sparse_available(), reason="SciPy not installed")
def test_sparse_alignment_is_as_accurate_as_dense():
    d, truth = _sim(taxa=5, root_len=150, seed=1, indel_rate=2.0, bmax=0.3)
    a = progressive.align(d["seqs"], Alphabet.DNA, sparse=False, guide_seed=0)
    b = progressive.align(d["seqs"], Alphabet.DNA, sparse=True, guide_seed=0)
    sa = evaluate.compare_to_reference(a, truth).sp
    sb = evaluate.compare_to_reference(b, truth).sp
    assert abs(sa - sb) < 0.05      # usually identical; never materially worse


def _sparse_choice(monkeypatch, seqs, **kw):
    """The ``sparse`` argument align() hands to the consistency transform."""
    seen = []
    real = progressive.consistency_transform

    def spy(*a, **k):
        seen.append(k.get("sparse"))
        return real(*a, **k)

    monkeypatch.setattr(progressive, "consistency_transform", spy)
    progressive.align(seqs, Alphabet.DNA, estimate=False, **kw)
    return seen


def test_align_uses_the_exact_transform_whenever_it_fits(monkeypatch):
    """The sparse transform is an approximation that costs accuracy on divergent
    sequences, so long sequences alone must not trigger it (0.5.10 switched to it
    above 350 residues whenever SciPy was installed)."""
    d, _truth = _sim(taxa=3, root_len=400, seed=2)
    assert _sparse_choice(monkeypatch, d["seqs"]) == [False]


@pytest.mark.skipif(not progressive.sparse_available(), reason="SciPy not installed")
def test_align_uses_the_sparse_transform_only_when_dense_does_not_fit(monkeypatch):
    d, _truth = _sim(taxa=3, root_len=400, seed=2)
    n, L = 3, max(len(s) for _, s in d["seqs"])
    budget = (progressive._posterior_gb(n, L, sparse=True) + progressive._posterior_gb(n, L)) / 2
    assert _sparse_choice(monkeypatch, d["seqs"], consistency_mem_gb=budget) == [True]


def test_memory_estimate_counts_both_generations():
    """The dense transform holds the new posteriors while the old are still
    referenced; v0.2's estimate counted only one set and so understated it."""
    one_pair = (500 ** 2) * 8 / 1e9
    assert progressive._posterior_gb(2, 500) == pytest.approx(2 * one_pair)
    assert progressive._posterior_gb(2, 500, sparse=True) < progressive._posterior_gb(2, 500)


# --------------------------------------------------------------------------- #
# Truth mode in the GUI
# --------------------------------------------------------------------------- #

def test_truth_mode_recomputes_when_the_alignment_changes():
    """The reference outlives a realignment — that is the whole teaching loop:
    align, look, realign differently, watch the correctness track change."""
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.engines import BuiltinProgressive
    from craic.gui.app import CraicWindow, _TRUTH_TRACK

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    d, truth = _sim(taxa=5, root_len=120, seed=5, indel_rate=3.0, bmax=0.4)

    width = max(len(s) for _, s in d["seqs"])
    win._set_alignment(Alignment([n for n, _ in d["seqs"]],
                                 [s.ljust(width, "-") for _, s in d["seqs"]],
                                 Alphabet.DNA))
    win.reference = truth
    aln = BuiltinProgressive().align(d["seqs"], Alphabet.DNA)
    win._set_alignment(aln)                      # "press Align"

    assert win.truth is not None
    assert win.truth.col_correct.size == aln.length
    win.track_combo.setCurrentIndex(_TRUTH_TRACK)
    assert win.canvas.score_label == "Reference correctness"
    assert len(win.canvas.nt_scores) == aln.length


def test_truth_mode_without_a_reference_explains_itself():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.gui.app import CraicWindow, _TRUTH_TRACK

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    _d, truth = _sim(taxa=4, root_len=80)
    win._set_alignment(truth)
    errors = []
    win._on_track_error = lambda msg: errors.append(msg)
    win.track_combo.setCurrentIndex(_TRUTH_TRACK)
    assert errors and "reference" in errors[0].lower()


def test_a_reference_for_other_data_is_dropped_quietly_on_reload():
    """Opening different sequences invalidates the reference; realigning does
    not. Only the first should clear it, and without a popup."""
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.gui.app import CraicWindow

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    _d, truth = _sim(taxa=4, root_len=80, seed=1)
    win._set_alignment(truth)
    win.reference = truth
    win._rescore_against_reference()
    assert win.truth is not None

    _d2, other = _sim(taxa=4, root_len=80, seed=42)
    win._set_alignment(other)                    # different sequences entirely
    assert win.reference is None and win.truth is None


def test_generated_dataset_loads_unaligned_and_keeps_the_answer():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication, QDialog

    from craic.gui import app as app_mod
    from craic.gui.app import CraicWindow

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()

    class FakeDialog:
        def __init__(self, parent=None):
            pass

        def exec(self):
            return QDialog.Accepted

        def values(self):
            return {"taxa": 5, "root_len": 100, "seed": 3, "indel_rate": 2.0,
                    "bmax": 0.3, "rate_alpha": 0.0}

    original = app_mod.SimulateDialog
    app_mod.SimulateDialog = FakeDialog
    try:
        win._generate_dataset()
    finally:
        app_mod.SimulateDialog = original

    assert win.aln is not None and win.aln.n_seqs == 5
    assert win.aln.meta.get("unaligned") is True
    assert win.reference is not None              # the answer is held, ready
    # …and scored at once, against the unaligned sequences. Holding it without
    # scoring left truth mode refusing to turn on, which reads as the program
    # having forgotten the dataset it had just generated.
    assert win.truth is not None
    assert win.doc.truth_cells is not None
    assert win._paint_truth_cells() is True


def test_a_cheap_track_survives_an_edit_and_the_numbers_move():
    """Truth mode has to stay live while you work: the whole point of holding the
    answer is watching an edit make the alignment better or worse."""
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic import editing
    from craic.engines import BuiltinProgressive
    from craic.gui.app import CraicWindow, _TRUTH_TRACK

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    d, truth = _sim(taxa=6, root_len=120, seed=5, indel_rate=3.0, bmax=0.4)
    win._set_alignment(BuiltinProgressive().align(d["seqs"], Alphabet.DNA))
    win._set_reference(truth, "simulated truth")
    assert win.track_combo.currentIndex() == _TRUTH_TRACK

    moved = False
    for col in range(5, 25):
        new = editing.move_residue(win.aln, 0, col=col, direction=1)
        if new is None:
            continue
        before = win.truth.sp
        win._apply_edit(new, "nudge")
        moved = True
        # the track is still showing, and re-scored against the reference
        assert win.track_combo.currentIndex() == _TRUTH_TRACK
        assert win.canvas.score_label == "Reference correctness"
        assert len(win.canvas.nt_scores) == win.aln.length
        assert win.truth is not None and win.truth.sp != before
        break
    assert moved, "no residue could be nudged in this alignment"
    assert "SP=" in win.sel_lbl.text()


def test_an_expensive_track_still_resets_on_an_edit():
    """Reliability takes seconds to minutes, so re-running it on every nudge
    would make editing unusable. It resets and is recomputed on request."""
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.gui.app import CraicWindow

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    _d, truth = _sim(taxa=4, root_len=80)
    win._set_alignment(truth)
    win.track_combo.setCurrentIndex(1)                     # Conservation, cheap
    win._set_alignment(truth)
    assert win.track_combo.currentIndex() == 1             # kept

    win.track_combo.blockSignals(True)
    win.track_combo.setCurrentIndex(2)                     # Reliability, expensive
    win.track_combo.blockSignals(False)
    win._set_alignment(truth)
    assert win.track_combo.currentIndex() == 0             # reset


def test_truth_track_is_dropped_if_the_reference_stops_applying():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.gui.app import CraicWindow, _TRUTH_TRACK

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    _d, truth = _sim(taxa=4, root_len=80, seed=1)
    win._set_alignment(truth)
    win._set_reference(truth, "truth")
    assert win.track_combo.currentIndex() == _TRUTH_TRACK

    _d2, other = _sim(taxa=4, root_len=80, seed=42)        # different sequences
    win._set_alignment(other)
    assert win.track_combo.currentIndex() == 0
    assert win.truth is None


# --------------------------------------------------------------------------- #
# Seeing the answer itself (v0.5): per-residue correctness, the column
# inspector, and the read-only view of the true alignment
# --------------------------------------------------------------------------- #

# One residue out of place, and nothing else. s3's G is put in column 2, where
# the other two have a C, instead of column 3 where it belongs.
_INF = ["ACGT", "ACGT", "AG-T"]
_TRUE = ["ACGT", "ACGT", "A-GT"]
_IDS = ["s1", "s2", "s3"]


def _tiny(rows):
    return Alignment(list(_IDS), list(rows), Alphabet.DNA)


def test_cell_correctness_names_the_residue_that_broke_the_column():
    cells = evaluate.cell_correct_fraction(_INF, _TRUE)
    # Column 1 (0-based): the two C's agree with each other but not with the
    # intruder, which agrees with nobody.
    assert cells[0][1] == pytest.approx(0.5)
    assert cells[1][1] == pytest.approx(0.5)
    assert cells[2][1] == pytest.approx(0.0)
    assert np.isnan(cells[2][2])                   # a gap scores nothing
    assert np.all(cells[:, 0] == 1.0)              # a correct column is all 1s


def test_a_columns_score_is_the_mean_of_its_cells():
    """The per-residue view must not be able to disagree with the track above it."""
    d, truth = _sim(taxa=5, root_len=150, seed=11, indel_rate=2.0, bmax=0.5)
    from craic.engines import BuiltinProgressive
    aln = BuiltinProgressive().align(d["seqs"], Alphabet.DNA)
    true_rows = evaluate.matched_rows(aln, truth)

    cells = evaluate.cell_correct_fraction(aln.rows, true_rows)
    cols = evaluate.col_correct_fraction(aln.rows, true_rows)
    for c in range(aln.length):
        column = cells[:, c]
        if np.all(np.isnan(column)):
            assert np.isnan(cols[c])
        else:
            assert np.nanmean(column) == pytest.approx(cols[c])


def test_the_inspector_states_the_correct_answer_for_a_wrong_column():
    e = evaluate.explain_column(_INF, _TRUE, 1, ids=_IDS)
    assert not e.correct
    assert e.true_col == 1
    assert [p.seq_id for p in e.agree] == ["s1", "s2"]
    intruder, = e.intruders
    assert intruder.seq_id == "s3" and intruder.char == "G"
    assert intruder.target_col == 2 and intruder.offset == 1      # one to the right
    assert not e.missing


def test_a_column_can_score_perfectly_and_still_be_incomplete():
    """Column 2 asserts only true pairs, so it scores 1.0 — but a residue that
    belongs in it was put elsewhere. Saying so is the point of the inspector."""
    assert evaluate.col_correct_fraction(_INF, _TRUE)[2] == pytest.approx(1.0)
    e = evaluate.explain_column(_INF, _TRUE, 2, ids=_IDS)
    assert not e.correct
    missing, = e.missing
    assert missing.seq_id == "s3" and missing.col == 1            # where it went
    assert not e.intruders


def test_a_correct_column_says_so():
    e = evaluate.explain_column(_INF, _TRUE, 0, ids=_IDS)
    assert e.correct and not e.intruders and not e.missing
    assert len(e.agree) == 3
    assert "exactly right" in e.describe()


def test_truth_colouring_and_the_inspector_follow_an_edit():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication
    from craic.gui.app import CraicWindow, _COLOUR_TRUTH

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    win._set_alignment(_tiny(_INF))
    win._set_reference(_tiny(_TRUE), "truth")
    win.color_combo.setCurrentIndex(_COLOUR_TRUTH)
    assert win.canvas.color_mode == "truth"
    assert win.doc.truth_cells[2][1] == pytest.approx(0.0)

    win._jump_to_column(1)
    assert win._inspect_col == 1

    # Fix it by hand: s3's G moves into the column where it belongs.
    win._set_alignment(_tiny(["ACGT", "ACGT", "A-GT"]))
    assert win.canvas.color_mode == "truth"          # the mode survives the edit
    assert np.nanmin(win.doc.truth_cells) == 1.0     # …and now everything is right
    e = win.doc.explain_column(1)
    assert e.correct


def test_truth_colouring_is_refused_without_an_answer():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication
    from craic.gui.app import CraicWindow, _COLOUR_RESIDUE, _COLOUR_TRUTH

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    win._set_alignment(_tiny(_INF))
    win.color_combo.blockSignals(True)
    win.color_combo.setCurrentIndex(_COLOUR_TRUTH)
    win.color_combo.blockSignals(False)
    win._paint_truth_cells()                        # the guard, without the dialog
    assert win.canvas.color_mode != "truth"
    win.color_combo.setCurrentIndex(_COLOUR_RESIDUE)
    assert win.canvas.color_mode == "residue"


def test_the_true_alignment_can_be_shown_and_is_read_only():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication
    from craic.gui.app import CraicWindow

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    win._set_alignment(_tiny(_INF))
    win._set_reference(_tiny(_TRUE), "truth")
    assert win.show_truth_act.isEnabled()

    win.show_truth_act.setChecked(True)
    assert win.showing_truth
    assert win.canvas.aln.rows == _TRUE            # the answer is on screen…
    assert win.doc.alignment.rows == _INF          # …but the document is untouched
    assert win.canvas.read_only
    assert not win.align_btn.isEnabled()

    win._set_dirty(False)
    win._apply_edit(_tiny(["A-CGT", "A-CGT", "AG--T"]), "refused")
    assert win.doc.alignment.rows == _INF
    assert not win.doc.dirty                       # a refused edit is not a change

    win.show_truth_act.setChecked(False)
    assert not win.showing_truth
    assert win.canvas.aln.rows == _INF
    assert not win.canvas.read_only


def test_forgetting_the_reference_takes_the_answer_off_the_screen():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication
    from craic.gui.app import CraicWindow, _COLOUR_RESIDUE, _COLOUR_TRUTH

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    win._set_alignment(_tiny(_INF))
    win._set_reference(_tiny(_TRUE), "truth")
    win.color_combo.setCurrentIndex(_COLOUR_TRUTH)
    win.show_truth_act.setChecked(True)

    win._clear_reference()
    assert not win.showing_truth
    assert win.canvas.color_mode == "residue"
    assert win.color_combo.currentIndex() == _COLOUR_RESIDUE
    assert not win.show_truth_act.isEnabled()
    assert win.doc.explain_column(1) is None


# --------------------------------------------------------------------------- #
# Genetic codes (v0.5.1): the codon views and "align as protein" are only as
# right as the translation table they use
# --------------------------------------------------------------------------- #

def test_the_code_list_covers_the_ones_people_actually_need():
    from craic.domain import genetic_codes, code_name
    codes = dict(genetic_codes())
    assert codes[1] == "Standard"
    assert "Vertebrate Mitochondrial" in codes[2]
    assert "Mitochondrial" in codes[4]          # also Mycoplasma / Spiroplasma
    assert 11 in codes                           # bacterial
    assert genetic_codes()[0][0] == 1            # the standard code is offered first
    assert code_name(2) == codes[2]
    assert code_name(999) == "table 999"         # unknown ids degrade, not raise


def test_tga_is_a_stop_or_a_tryptophan_depending_on_the_table():
    from craic.domain import translate_codon
    assert translate_codon("TGA", 1) == "*"
    assert translate_codon("TGA", 4) == "W"      # Mycoplasma / Spiroplasma
    assert translate_codon("TGA", 2) == "W"      # vertebrate mitochondrial
    assert translate_codon("ATA", 1) == "I"
    assert translate_codon("ATA", 2) == "M"


def _coding(rows, table):
    from craic.domain import CodingSpec
    return Alignment(["s%d" % i for i in range(len(rows))], list(rows),
                     Alphabet.DNA, coding=CodingSpec(frame=0, table=table))


def test_internal_stops_are_found_under_the_wrong_code_and_not_the_right_one():
    """A Mollicute gene read with the standard code fills with stop codons."""
    rows = ["ATGAAATGGTGAAAATGAGGTTTATAA"] * 2      # TGA twice, TAA terminal
    assert len(_coding(rows, 1).internal_stops()) == 4     # 2 per sequence
    assert _coding(rows, 4).internal_stops() == []


def test_a_terminal_stop_is_not_an_internal_stop():
    assert _coding(["ATGAAATAA"], 1).internal_stops() == []


def test_trailing_gaps_do_not_make_a_short_sequence_stop_early():
    """A sequence padded out to the alignment width still ends at its own last
    codon, so its real terminal stop is not reported."""
    rows = ["ATGAAATAA---", "ATGAAAAAATAA"]
    assert _coding(rows, 1).internal_stops() == []


def test_the_codon_engine_honours_the_table():
    from craic.engines import BuiltinProgressive, CodonAware
    recs = [("a", "ATGTGGTGAAAA"), ("b", "ATGTGGTGAAAA")]
    out1 = CodonAware(BuiltinProgressive(), table=1).align(recs, Alphabet.DNA)
    out4 = CodonAware(BuiltinProgressive(), table=4).align(recs, Alphabet.DNA)
    assert out1.coding.table == 1 and out4.coding.table == 4
    assert "*" in out1.translated_rows()[0]        # TGA read as a stop
    assert "*" not in out4.translated_rows()[0]    # TGA read as tryptophan


def test_a_stop_codon_carries_no_homology_signal():
    """Why the table matters: '*' is not one of the twenty residues, so it
    encodes to the wildcard and contributes nothing to the posterior."""
    from craic.accel import emission_model
    m = emission_model("protein")
    assert m.encode("*") == m.encode("X") == [m.index["X"]]


def test_cli_exposes_the_codes_and_passes_the_table_through(tmp_path, capsys):
    src = tmp_path / "in.fasta"
    src.write_text(">a\nATGTGGTGAAAA\n>b\nATGTGGTGAAAA\n")
    out = tmp_path / "out.fasta"
    from craic.cli import main
    with pytest.raises(SystemExit):
        main(["codes"])
    listed = capsys.readouterr().out
    assert "Standard" in listed and "Mitochondrial" in listed

    with pytest.raises(SystemExit):
        main(["align", str(src), "--codon", "--code", "4", "-o", str(out)])
    from craic import io
    assert io.load_alignment(str(out)).n_seqs == 2


def test_gui_offers_the_code_and_applies_it():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication, QMessageBox
    from craic.gui.app import CraicWindow

    app = QApplication.instance() or QApplication([])      # noqa: F841
    seen = []
    QMessageBox.warning = staticmethod(lambda *a, **k: seen.append(a) or QMessageBox.Ok)

    win = CraicWindow()
    win._set_alignment(Alignment(["a", "b"], ["ATGAAATGGTGAAAATAA"] * 2, Alphabet.DNA))
    assert win.code_combo.count() > 10

    win.coding_chk.setChecked(True)
    assert win.aln.coding.table == 1
    assert seen, "an internal stop under the standard code should be reported"

    win.code_combo.setCurrentIndex(win.code_combo.findData(4))
    assert win.aln.coding.table == 4
    assert win.aln.internal_stops() == []
    # TGA is now tryptophan; the trailing TAA is a real stop under both codes
    assert win.aln.translated_rows()[0] == "MKWWK*"


def test_the_default_engine_is_the_builtin_one():
    """The start-up engine is the built-in one, on purpose.

    It is the only engine guaranteed to be installed, so it is the only one that
    can be the default on every machine, and it is the engine the reliability
    overlay, posterior explorer, sandbox and perturbation ensemble all read their
    posteriors from. On BAliBASE it is level with MAFFT and ahead of Clustal
    Omega, behind MUSCLE and ProbCons.
    """
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.gui.app import CraicWindow

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    assert win.engine_combo.currentData().key == "builtin"


def test_pinned_columns_can_be_set_survive_an_edit_and_drive_the_realign():
    """Pinning was documented and menu-named but never implemented: nothing in
    the package wrote ``meta['anchors']``, so the realign always fell through to
    a hidden fallback on the current selection."""
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.domain import Alignment
    from craic.gui.app import CraicWindow

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    aln = Alignment(["a", "b", "c"],
                    ["ACGTACGTAC", "ACGTACGTAC", "ACGTTCGTAC"], Alphabet.DNA)
    win._set_alignment(aln)

    win._set_pins([2, 3, 7])
    assert win._pins() == {2, 3, 7}
    assert win.aln.meta["anchors"] == [2, 3, 7]        # JSON-safe, so sessions keep it
    assert win.canvas.pins == {2, 3, 7}

    # out-of-range pins are refused rather than stored
    win._set_pins([1, 99])
    assert win._pins() == {1}

    # …and they survive an edit that keeps the columns
    win._set_pins([2, 3, 7])
    win._commit(Alignment(list(aln.ids), list(aln.rows), Alphabet.DNA), "no-op edit")
    assert win._pins() == {2, 3, 7}

    # …but an edit that changes the width drops them, because they have moved and
    # a pin pointing at whatever slid into its place is worse than no pin
    win._commit(Alignment(["a", "b", "c"], ["ACGTA", "ACGTA", "ACGTT"], Alphabet.DNA),
                "trimmed")
    assert win._pins() == set()

    win._clear_pins()
    assert win._pins() == set() and "anchors" not in win.aln.meta


def test_an_anchored_realign_moves_the_pins_with_their_columns():
    """The blocks between anchors change width when they are re-solved, so the
    anchors land at new indices. Keeping the old ones would silently re-pin
    whatever moved into that position."""
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication, QMessageBox

    from craic.engines import BuiltinProgressive
    from craic.gui.app import CraicWindow

    app = QApplication.instance() or QApplication([])      # noqa: F841
    QMessageBox.information = staticmethod(lambda *a, **k: QMessageBox.Ok)
    win = CraicWindow()
    d, _truth = _sim(taxa=5, root_len=120, seed=4, indel_rate=2.0, bmax=0.4)
    win._set_alignment(BuiltinProgressive().align(d["seqs"], Alphabet.DNA))

    pins = [5, 6, 7, 20, 21]
    win._set_pins(pins)
    before = [[r[c] for r in win.aln.rows] for c in pins]
    win._anchored_realign()

    after_cols = sorted(win._pins())
    assert len(after_cols) == len(pins)
    # the pinned columns still hold exactly the residues they held before
    assert [[r[c] for r in win.aln.rows] for c in after_cols] == before


def test_a_protein_only_engine_is_refused_before_it_runs_not_after():
    """ProbCons cannot align nucleotides. That must be said by the interface,
    not by a traceback escaping a worker thread."""
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication, QMessageBox

    from craic.engines import ProbConsEngine
    from craic.gui.app import CraicWindow

    app = QApplication.instance() or QApplication([])      # noqa: F841
    said = []
    QMessageBox.information = staticmethod(lambda *a, **k: (said.append(a[2]), QMessageBox.Ok)[1])
    win = CraicWindow()
    d, _truth = _sim(taxa=4, root_len=90, seed=2)
    from craic.engines import BuiltinProgressive
    win._set_alignment(BuiltinProgressive().align(d["seqs"], Alphabet.DNA))

    pc = ProbConsEngine()
    assert pc.supports(Alphabet.PROTEIN) and not pc.supports(Alphabet.DNA)

    # the combo refuses to sit on an engine that cannot take the data
    idx = next((i for i in range(win.engine_combo.count())
                if getattr(win.engine_combo.itemData(i), "key", "") == "probcons"), None)
    if idx is not None:
        win.engine_combo.setCurrentIndex(idx)
        win._sync_engine_combo()
        assert win.engine_combo.currentData().key != "probcons"

    # and asking directly explains rather than raising
    win.engine_combo.addItem(pc.label, pc)
    win.engine_combo.setCurrentIndex(win.engine_combo.count() - 1)
    win._align()
    assert said and "protein sequences only" in said[-1]
    assert win._align_token is None                 # nothing was ever started


def test_align_as_protein_unlocks_a_protein_only_engine():
    """Ticking “align as protein” translates the nucleotides before the engine
    sees them, so ProbCons is perfectly usable on coding DNA — it is handed amino
    acids. Judging the engine by the document's own alphabet ruled out exactly
    the case the checkbox exists to enable."""
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from craic.domain import CodingSpec
    from craic.engines import BuiltinProgressive, ProbConsEngine
    from craic.gui.app import CraicWindow

    app = QApplication.instance() or QApplication([])      # noqa: F841
    win = CraicWindow()
    d, _truth = _sim(taxa=4, root_len=120, seed=3)
    aln = BuiltinProgressive().align(d["seqs"], Alphabet.DNA)
    aln.coding = CodingSpec(frame=0, table=1)
    win._set_alignment(aln)
    win.engine_combo.addItem(ProbConsEngine().label, ProbConsEngine())
    idx = win.engine_combo.count() - 1
    model = win.engine_combo.model()

    win.codon_chk.setChecked(False)
    win._sync_engine_combo()
    assert win._engine_alphabet() == Alphabet.DNA
    assert not model.item(idx).isEnabled()

    win.codon_chk.setChecked(True)
    win._sync_engine_combo()
    assert win._engine_alphabet() == Alphabet.PROTEIN
    assert model.item(idx).isEnabled()


def test_engine_parameters_refuse_values_the_aligner_would_reject():
    """ProbCons takes 0–1000 refinement passes. Typing 10000 used to be accepted
    by the dialog and rejected by ProbCons after the run had started."""
    pytest.importorskip("PySide6")
    from PySide6.QtGui import QValidator
    from PySide6.QtWidgets import QApplication

    from craic.engines import ProbConsEngine
    from craic.gui.dialogs import EngineParamsDialog

    app = QApplication.instance() or QApplication([])      # noqa: F841
    dlg = EngineParamsDialog(ProbConsEngine(), {})
    p, w = dlg._widgets["iterations"]
    assert (p.lo, p.hi) == (0, 1000)

    v = w.validator()
    assert v.validate("10000", 5)[0] == QValidator.Invalid
    assert v.validate("1000", 4)[0] == QValidator.Acceptable
    assert v.validate("-5", 2)[0] == QValidator.Invalid

    # paste bypasses the validator on some platforms, so the read clamps too
    w.setText("99999")
    assert dlg.values()["iterations"] == 1000
