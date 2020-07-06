call cd src

nmake -f Makefile.win cleanall

REM Generate Python3.6
call ..\env36\Scripts\activate.bat
call nmake -f Makefile.win build_ext
call python setup.py bdist_wheel
call deactivate

call ..\env37\Scripts\activate
call nmake -f Makefile.win build_ext
call python setup.py bdist_wheel
call deactivate

call ..\env38\Scripts\activate
call nmake -f Makefile.win build_ext
call python setup.py bdist_wheel
call deactivate

call cd ..
