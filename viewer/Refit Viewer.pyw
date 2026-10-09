# Double-click to open the refit viewer in its own window, without a console.
import os, runpy, sys

sys.argv = [os.path.join(os.path.dirname(os.path.abspath(__file__)), "server.py")] + sys.argv[1:]
runpy.run_path(sys.argv[0], run_name="__main__")
