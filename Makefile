install:
	python3 -m pip install --upgrade pip
	python3 -m pip install --upgrade setuptools	
	python3 -m pip install --upgrade build
	python3 -m pip install --upgrade colorama
	python3 -m build --no-isolation .
	python3 -m pip install dist/krisp-1.0.0-py3-none-any.whl --force-reinstall

uninstall:
	python3 -m pip uninstall krisp

