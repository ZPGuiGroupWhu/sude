import pathlib
import unittest


try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11
    import tomli as tomllib


PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]


class TestPackagingMetadata(unittest.TestCase):
    def test_default_dependencies_include_numba(self):
        pyproject = tomllib.loads((PROJECT_ROOT / "pyproject.toml").read_text())

        dependencies = pyproject["project"]["dependencies"]

        self.assertTrue(
            any(dependency.startswith("numba") for dependency in dependencies),
            "numba should be installed by default with pip install sude",
        )

    def test_readme_documents_default_acceleration(self):
        readme = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")

        self.assertIn("Numba-accelerated kernels are installed by default", readme)
        self.assertNotIn("sude[accelerate]", readme)


if __name__ == "__main__":
    unittest.main()
