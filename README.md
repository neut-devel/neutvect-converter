# neutvect-converter
Tools for converting neutvect ROOT files to another format (Currently just something approximating [NuHepMC](https://github.com/NuHepMC/Spec) version 0.9).

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

## Usage

```
$ neutvect-converter -?
[USAGE]: neutvect-converter
  -i <neutvect.root>       : neutvect file to read
  -N <NMax>                : Process at most <NMax> events
  -o <neut.hepmc3>         : hepmc3 file to write
  -f <flux_file,flux_hist> : ROOT flux histogram to use to
  -z                       : Write to .gz compress ASCII file
  -G                       : -f argument should be interpreted as being in GeV
```

For the majority of files -f and -G options are not required as the input neutvect file will contain enough information to calculate the flux-averaged total cross section, but if you really need to pass a flux, you can.
