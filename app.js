import { cardan2R3D, readXYZFile, parseXYZData } from "./lib/utils.js";

readXYZFile("./data/trajWGS84.xyz").then((d) =>
  console.log(parseXYZData(d))
);
