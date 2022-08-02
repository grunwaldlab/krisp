install:
	python3 -m pip install --upgrade pip
	python3 -m pip install --upgrade setuptools	
	python3 -m pip install --upgrade build
	python3 -m pip install --upgrade colorama
	python3 -m build .
	python3 -m pip install --editable .

uninstall:
	python3 -m pip uninstall krisp

test:
	python3 -m unittest discover
