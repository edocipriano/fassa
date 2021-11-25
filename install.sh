rm -rvf build dist fassa.egg-info
python3 setup.py bdist_wheel
pip3 install dist/*
