"""Startup dependency check: installs what is missing with pip.

ensure() is called by refit.py and viewer/server.py before their first
third-party import. It only touches the network when something is missing.
Set REFIT_NO_AUTOINSTALL=1 to turn the installation off (the check still runs).
"""
import importlib.util
import os
import subprocess
import sys

# module name -> pip package, required
REQUIRED = {"numpy": "numpy", "scipy": "scipy", "OCP": "cadquery-ocp"}
# optional: mesh repair (refit.py) and the viewer window
OPTIONAL = {"pymeshlab": "pymeshlab", "webview": "pywebview"}


def _missing(table):
    return [pkg for mod, pkg in table.items() if importlib.util.find_spec(mod) is None]


def _pip(pkgs):
    exe = sys.executable
    if os.path.basename(exe).lower() == "pythonw.exe":  # .pyw launch: no console
        exe = os.path.join(os.path.dirname(exe), "python.exe")
    flags = 0x08000000 if os.name == "nt" else 0  # CREATE_NO_WINDOW
    r = subprocess.run([exe, "-m", "pip", "install", *pkgs], creationflags=flags)
    return r.returncode == 0


def ensure(optional=True):
    """Install missing packages. Exits if a required one cannot be installed."""
    req, opt = _missing(REQUIRED), (_missing(OPTIONAL) if optional else [])
    if not (req or opt):
        return
    if os.environ.get("REFIT_NO_AUTOINSTALL"):
        if req:
            sys.exit("Missing dependencies: pip install " + " ".join(req))
        return
    for group, name in ((req, "required"), (opt, "optional")):
        if not group:
            continue
        print(f"[deps] installing {name} packages: {' '.join(group)}", flush=True)
        if not _pip(group):
            if name == "required":
                sys.exit("Could not install: pip install " + " ".join(group))
            print("[deps] optional install failed, continuing without", flush=True)
    importlib.invalidate_caches()
    if _missing(REQUIRED):
        sys.exit("Missing dependencies after install: " + " ".join(_missing(REQUIRED)))
    # OCP can install yet fail to load without VTK (see README)
    if importlib.util.find_spec("OCP") is not None:
        try:
            import OCP  # noqa: F401
        except ImportError as ex:
            if "DLL" in str(ex) and _pip(["vtk==9.6.2"]):
                importlib.invalidate_caches()


if __name__ == "__main__":
    ensure()
    print("[deps] ok")
