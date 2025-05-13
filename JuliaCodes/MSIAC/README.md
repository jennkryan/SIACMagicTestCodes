## Getting started

Install MSIAC as a Julia package in our local distribution

Go to terminal and run 
```bash
$ cd path_to_this_folder
$ make dev
```

### If you got an error, it might be that you are not in shell. Make sure  that 
1. julia is a command to launch Julia
2. You are running on bash 

## Unit Testing 
```bash
$ make test
```

It should pass a bunch of  tests


***If everything worked, now your local Julia has a MSIAC package in its package list. 
You should be able to run a Julia session from wherever in your machine and load MSIAC***
```julia
using MSIAC  (or import MSIAC)
```

## Resources 
1. In the [test folder](https://gitlab.com/msiac-tool/MSIAC/-/tree/main/test), there are several examples that you can checkout
2. The [tutorials folder](https://gitlab.com/msiac-tool/MSIAC/-/tree/main/tutorials) also has a collection of examples using different filter

