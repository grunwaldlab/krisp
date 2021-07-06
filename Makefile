install:
	python3 -m pip install  --upgrade --user pip
	python3 -m pip install --upgrade --user build
	python3 -m build --no-isolation .
	python3 -m pip install dist/krisp-1.0.0-py3-none-any.whl

uninstall:
	python3 -m pip uninstall krisp

