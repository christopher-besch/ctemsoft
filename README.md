# Chris' Attempt
- `sudo docker build -t chrisbesch/ctemsoft . && sudo docker run -v ./out:/python_src/out chrisbesch/ctemsoft`
- https://ctem.web.cmu.edu
- https://ctem.web.cmu.edu/frames.html
- `docker run --net host --name lea-fortran -ti -v .:/code debian`
```
    1  ls
    2  set -o vi
    3  ls
    4  cd /code
    5  ls
    6  apt update
    7  ls
    8  cd src/
    9  ls
   10  make
   11  apt install make
   12  make
   13  make
   14  f90
   15  apt install build-essential
   16  apt install gfortran
   17  make
   18  apt search lapack
   19  apt install liblapack3
   20  apt search blas
   21  apt install libblas3 liblapack3 liblapack-dev libblas-dev
   22  apt insatll libblas3
   23  apt install libblas3
   24  which lapack-3.0
   25  ls /usr/local/lib
   26  ls /usr/lib
   27  dpkg -L liblapack3
   28  cd /usr/lib/x86_64-linux-gnu/lapack
   29  ls
   30  file *
   31  apt insatll file
   32  apt install file
   33  file *
   34  ls
   35  dpkg -L libblas3
   36  dpkg -L liblapack3
   37  apt install libfftw3
   38  apt install libfftw3-3
   39  apt search fftw
   40  apt search libfftw3
   41  apt install libfftw3
   42  apt install libfftw3-dev
   43  apt install libfftw3-bin libfftw3-double3
   44  dpkg -L libfftw3-dev
   45  dpkg -L libfftw3-dev | grep libsfftw.a
   46  dpkg -L libfftw3-double3 | grep libsfftw.a
   47  dpkg -L libfftw3-single3 | grep libsfftw.a
   48  dpkg -L libfftw3-single3
   49  cd
   50  cd code
   51  cd /code/src
   52  ls
   53  make
   54  which f90
   55  f90
   56  f77
   57  gfortran --version
   58  gfortran --help
   59  make
   60  make
   61  make
   62  make
   63  make
   64  make
   65  make
   66  make
   67  make CHAP3
   68  make chap3
   69  cd /usr/bin/ld: cannot find /usr/lib/x86_64-linux-gnu/lapack/
   70  ls
   71  cd /usr/lib/x86_64-linux-gnu/lapack/
   72  ls
   73  make chap3
   74  cd
   75  cd /code/src
   76  make chap3
   77  which libtem.a
   78  make chap3
   79  ls
   80  make chap3
   81  make LIBTEM
   82  make libtem.a
   83  make ctf.f90
   84  make libtem.a
   85  clear
   86  make libtem.a
   87  make libtem.a
   88  make libtem.a
   89  make libtem.a
   90  make libtem.a
   91  make libtem.a
   92  apt install gfortran-doc
   93  apt install gfortran-doc
   94  man gfortran
   95  apt install man
   96  man gfortran
   97  g90
   98  g95
   99  make libtem.a
  100  make libtem.a
  101  make clean
  102  make libtem.a
  103  clear
  104  make clean
  105  clear
  106  make libtem.a
  107  make clean
  108  clear
  109  make libtem.a
  110  make clean
  111  clear
  112  make libtem.a
  113  make clean
  114  clear
  115  make libtem.a
  116  make clean
  117  make libtem.a
  118  make clean
  119  make libtem.a
  120  clear
  121  make libtem.a
  122  clear
  123  make libtem.a
  124  clear
  125  make libtem.a
  126  clear
  127  make clean
  128  make libtem.a
  129  clear
  130  make clean
  131  clear
  132  make libtem.a
  133  gfortran --version
  134  make clean
  135  clear
  136  make libtem.a
  137  make libtem.a
  138  make clean
  139  make libtem.a
  140  make libtem.a
  141  make clean
  142  make libtem.a
  143  make clean
  144  make libtem.a
  145  make clean
  146  make libtem.a
  147  make clean
  148  make clean
  149  make libtem.a
  150  make clean
  151  make libtem.a
  152  make chap3
  153  echo $?
  154  make clean
  155  make libtem.a
  156  echo $?
  157  make chap3
  158  make all
  159  make chap3
  160  history | grep dpkg
  161  dpkg -L liblapack3
  162  find / -name \*.a -exec bash -c "nm --defined-only {} 2>/dev/null | grep 'cgetrf_' && echo {}" \;
  163  make chap3
  164  find / -name \*.a -exec bash -c "nm --defined-only {} 2>/dev/null | grep 'zdotc_' && echo {}" \;
  165  make clean
  166  make chap3
  167  make local.mod
  168  make
  169  echo $?
  170  ls
  171  ..
  172  ls
  173  cd ..
  174  ls
  175  ls exe/
  176  cd src/
  177  ls
  178  make all
  179  make clean
  180  ls
  181  make libtem.a
  182  make chap3
  183  find / -name \*.a -exec bash -c "nm --defined-only {} 2>/dev/null | grep 'zaxpy_' && echo {}" \;
  184  make chap3 | less
  185  apt install less
  186  make chap3 | less
  187  make chap3 2>&1 | less
  188  find / -name \*.a -exec bash -c "nm --defined-only {} 2>/dev/null | grep 'cgetri_' && echo {}" \;
  189  file /usr/lib/x86_64-linux-gnu/liblapack.a
  190  file /etc/alternatives/liblapack.a-x86_64-linux-gnu
  191  file /usr/lib/x86_64-linux-gnu/lapack/liblapack.a
  192  make chap3
  193  find / -name \*.a -exec bash -c "nm --defined-only {} 2>/dev/null | grep 'zscal_' && echo {}" \;
  194  make chap3
  195  echo $?
  196  ls
  197  ls ../exe
  198  cd ..
  199  ls
  200  cd exe
  201  ls
  202  file holz
  203  chmod +x *
  204  ls
  205  ll
  206  ls -l
  207  ./holz
  208  ./holz
  209  ls
  210  ./lens
  211  apt install gv
  212  ls
  213  ls
  214  ./lens
  215  ./lens
  216  ..
  217  cd ..
  218  cd src/
  219  cd ..
  220  cd exe/
  221  ls
  222  rm lens.ps
  223  ./lens
  224  ./lens
  225  ./lens
  226  ./lens
  227  ./lens
  228  ls
  229  ..
  230  cd ..
  231  cd src
  232  make lens
  233  make chap3
  234  cd ..
  235  cd exe
  236  ls
  237  ./lens
  238  cd ..
  239  cd src
  240  make chap3
  241  alias ..='cd ..'
  242  ..
  243  cd exe/
  244  ./lens
  245  ..
  246  cd src/
  247  make chap3
  248  ..
  249  cd exe/
  250  ./lens
  251  ./lens
  252  ..
  253  cd src
  254  make chap3
  255  ..
  256  cd exe/
  257  ./lens
  258  ..
  259  cd src/
  260  make chap3
  261  ..
  262  cd exe/
  263  ./lens
  264  echo $?
  265  ls
  266  ..
  267  cd src/
  268  make chap3
  269  ..c
  270  ..
  271  cd exe/
  272  ./lens
  273  ls
  274  ./holz
  275  ..
  276  cd src/
  277  make clean
  278  make libtem.a
  279  make chap3
  280  ..
  281  cd exe/
  282  ls
  283  ls -l
  284  apt install gdb
  285  ls
  286  gdb lens
  287  history
```
