# neutvect-converter
Tools for converting neutvect ROOT files to another format.

## build requirements

ROOT 6+
NEUT 5.5+
C++17 Compiler

## build like

```bash
cd /path/to/neutvect-converter
mkdir build; cd build
cmake ..
make install
```

## run like

```
/path/to/neutvect-converter/build/Linux/bin/neutvect-converter neutvect.root myout.hepmc3
```
