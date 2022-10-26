# helmert.js
JavaScript implementation of the Helmert transformation algorithm developed in Python by Sebastien Guillaume and Samuel Schwyn from INSIT at HEIG-VD. The goal of this project is to eventually integrate it in Biodivar, Julien Mercier doctoral thesis project whithin the Media Engineering institute, also at HEIG-VD.

## Warning
**As of right now, this project is neither functionnal nor architecturally sound, as we are still looking for a suitable SVD algorithm.** Right now, this stuff is built on Node, but the aim is to eventually have it running in the browser. 

## Get it running
To run all of the scripts, you will need:
* node.js
* Python 3

```bash
git clone https://github.com/MediaComem/helmert.js.git && cd helmert.js
npm install 
pip install matplotlib
```

## Scripts
`npm run dev`: Runs a demo of the JavaScript implementation using the files `traj.xyz` files in `data/`

`npm run py`: Runs a demo of the Python implementation using those same files.

## Test Protocol
The main of this project thus far is to try to replicate the results from the python program in the javascript implementation. To do so, we have implemented a test which does the following:
* Pseudo-ramdomly generates local and global datapoints in Python and outputs them in .xyz files in `data/generated_data/testLocalData.xyz` and `data/generated_data/testLocalData.xyz`
* The files are then fed to the Python algo and the resulting transformation is output to another file `data/generated_data/result.xyz`
* A Jest test compares the output of the JavaScript algorithm to the results in that file.

You can run the **as of now failing** test with the following command:
```bash
npm test
```
