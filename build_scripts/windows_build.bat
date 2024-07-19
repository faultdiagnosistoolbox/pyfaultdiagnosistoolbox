call cd src

nmake -f Makefile.win cleanall

call ..\..\env37\Scripts\activate
call nmake -f Makefile.win build_ext
call python setup.py bdist_wheel
call deactivate

call ..\..\env38\Scripts\activate
call nmake -f Makefile.win build_ext
call python setup.py bdist_wheel
call deactivate

call ..\..\env39\Scripts\activate
call nmake -f Makefile.win build_ext
call python setup.py bdist_wheel
call deactivate

call cd ..
