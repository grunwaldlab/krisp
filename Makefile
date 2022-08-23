install:
	python -m pip install --upgrade pip
	python -m pip install --upgrade setuptools
	python -m pip install --upgrade build
	python -m pip install --upgrade colorama
	python -m build .
	python -m pip install --editable .

uninstall:
	python -m pip uninstall krisp

test:
	python -m unittest discover
