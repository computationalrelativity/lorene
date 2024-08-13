# LORENE Installation on various systems in used by the CoRe collaboration

**WARNING** This folder is not part of the LORENE library but belongs to the `installation` branch of the [repo](https://github.com/computationalrelativity/lorene/wiki)

Installation notes for LORENE can be found [here](https://lorene.obspm.fr/install.html). In short, if you have cloned the library into a `lorene` folder, you do

```
cd lorene
export HOME_LORENE=$PWD # add this line to your bashrc, for later use
make # this uses a local_setting file in the current folder
make test
```

the important things are the `local_setting` file and the modules.

This folder contains examples of `local_setting` files and the list of modules for various systems used by CoRe collaboration. Please keep these files updated and use the following names

```
  <system>_<year>_<spec>.local_setting
  <system>_<year>_<spec>.modules
```

where

 * `<system>_<year>` indicates the machine and the year of a successful compilation.
 * `<spec>` is an optional argument in case a special setup is used. Comments about the latter should be included in the `modules` file.
