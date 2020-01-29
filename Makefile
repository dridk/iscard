SHELL=/bin/bash

test:
	py.test -vs --cov=iscard tests --cov-report term

black:
	black iscard


check: black style test

docs:
	make -f docs/Makefile html

# release:
# 	rm dist/*
# 	python setup.py sdist bdist_wheel
# 	twine upload dist/*
