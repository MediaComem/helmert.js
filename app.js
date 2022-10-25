import { Helmert } from "./lib/helmert.js";
import { eye, cardan2R3D, readXYZFile, parseXYZData } from "./lib/utils.js";

const helm = new Helmert();
helm
  .importFiles("./data/trajLocal.xyz", "./data/trajWGS84.xyz")
  .then(d => helm.estimateHelmertMinimum());

// readXYZFile("./data/trajWGS84.xyz").then((d) =>
//   console.log(parseXYZData(d))
// );
