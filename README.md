# jsondotrulo

Converts the .json output by Yosys to DOT-format and other graph descriptions.

### Init

`git submodule update --init --recursive`

### Build

* have cmake

```
mkdir build && cd build
cmake ../
make -j $(nproc)
```

### Running

```
build/jsondotrulo yosys_output.json --m_file=/path/to/m_prefix --dot_file=/path/to/dot_prefix
```

### Using graphviz

Do something like

```
neato -Tps -Goverlap="scale" yosys_output.gv -o yosys_output.ps
okular yosys_output.ps &
```

### Tests

hahaha
