import { Helmert } from "./src/helmert.js";

const helm = new Helmert();
helm.importFiles("./data/trajLocal.xyz", "./data/trajWGS84.xyz").then((d) => {
  console.log(helm.estimateHelmertMinimum().globalToLocal())
});
