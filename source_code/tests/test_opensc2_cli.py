from pathlib import Path

import pytest

import opensc2


@pytest.mark.parametrize("mode", ("recovery", "continuation"))
def test_parser_accepts_checkpoint_modes(mode):
    args = opensc2.parse_command_line_arguments(
        [
            "--no-head",
            "--io-path",
            "io.yaml",
            "--checkpoint",
            "state.h5",
            "--restart-mode",
            mode,
        ]
    )

    assert args.no_head is True
    assert args.io_path == "io.yaml"
    assert args.checkpoint == "state.h5"
    assert args.restart_mode == mode


def test_parser_preserves_normal_headless_mode():
    args = opensc2.parse_command_line_arguments(
        ["--no-head", "--io-path", "io.yaml"]
    )

    assert args.no_head is True
    assert args.checkpoint is None
    assert args.restart_mode is None


@pytest.mark.parametrize(
    "arguments",
    (
        ["--no-head", "--checkpoint", "state.h5"],
        ["--no-head", "--restart-mode", "recovery"],
        [
            "--checkpoint",
            "state.h5",
            "--restart-mode",
            "recovery",
        ],
    ),
)
def test_parser_rejects_incomplete_or_gui_checkpoint_requests(arguments):
    with pytest.raises(SystemExit) as error:
        opensc2.parse_command_line_arguments(arguments)

    assert error.value.code == 2


class FakeSimulation:
    def __init__(self, input_dir, events):
        self.input_dir = input_dir
        self.events = events
        self.dict_path = {}
        self.flag_start = False
        self.transient_input = {"SIMULATION": "cli_test"}

    def conductor_instance(self):
        self.events.append("conductor_instance")

    def simulation_folders_manager(self):
        self.events.append("simulation_folders_manager")

    def save_input_files(self):
        self.events.append("save_input_files")

    def conductor_initialization(self, gui):
        assert gui.simulation is self
        self.events.append("conductor_initialization")

    def conductor_solution(self, gui):
        assert gui.simulation is self
        self.events.append("conductor_solution")

    def conductor_post_processing(self):
        self.events.append("conductor_post_processing")


def configure_fake_runtime(monkeypatch, tmp_path):
    events = []
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()

    simulation_holder = {}

    def build_simulation(path):
        simulation = FakeSimulation(path, events)
        simulation_holder["simulation"] = simulation
        return simulation

    monkeypatch.setattr(
        opensc2,
        "load_io_path_from_yaml",
        lambda io_path: (input_dir, output_dir),
    )
    monkeypatch.setattr(opensc2, "Simulation", build_simulation)
    monkeypatch.setattr(opensc2.time, "time", lambda: 0.0)

    return events, simulation_holder, input_dir, output_dir


def test_normal_headless_run_does_not_restore_checkpoint(monkeypatch, tmp_path):
    events, holder, input_dir, output_dir = configure_fake_runtime(
        monkeypatch,
        tmp_path,
    )
    monkeypatch.setattr(
        opensc2,
        "read_checkpoint",
        lambda path: pytest.fail("Normal run must not read a checkpoint."),
    )
    monkeypatch.setattr(
        opensc2,
        "apply_checkpoint_to_runtime",
        lambda *args, **kwargs: pytest.fail(
            "Normal run must not restore a checkpoint."
        ),
    )

    opensc2.run_headless_simulation("io.yaml")

    simulation = holder["simulation"]
    assert simulation.input_dir == str(input_dir)
    assert simulation.dict_path["Main_dir"] == str(output_dir)
    assert simulation.flag_start is True
    assert events == [
        "conductor_instance",
        "simulation_folders_manager",
        "save_input_files",
        "conductor_initialization",
        "conductor_solution",
        "conductor_post_processing",
    ]


@pytest.mark.parametrize("mode", ("recovery", "continuation"))
def test_checkpoint_is_applied_after_initialization_and_before_solution(
    monkeypatch,
    tmp_path,
    mode,
):
    expected_mode = mode
    events, holder, _, _ = configure_fake_runtime(monkeypatch, tmp_path)
    checkpoint_path = tmp_path / "state.h5"
    checkpoint_path.touch()
    checkpoint = object()

    def read(path):
        assert Path(path) == checkpoint_path
        events.append("read_checkpoint")
        return checkpoint

    def apply(restored_checkpoint, simulation, *, mode):
        assert restored_checkpoint is checkpoint
        assert simulation is holder["simulation"]
        assert mode == expected_mode
        events.append(f"apply_checkpoint:{mode}")

    monkeypatch.setattr(opensc2, "read_checkpoint", read)
    monkeypatch.setattr(opensc2, "apply_checkpoint_to_runtime", apply)

    opensc2.run_headless_simulation(
        "io.yaml",
        checkpoint_path=checkpoint_path,
        restart_mode=mode,
    )

    assert events == [
        "read_checkpoint",
        "conductor_instance",
        "simulation_folders_manager",
        "save_input_files",
        "conductor_initialization",
        f"apply_checkpoint:{mode}",
        "conductor_solution",
        "conductor_post_processing",
    ]
