import { Helmert } from "./src/helmert.js";

const helm = new Helmert();
helm.importFiles("./data/trajLocal.xyz", "./data/trajWGS84.xyz").then((d) => {
  const localPoints = helm.estimateHelmertMinimum().globalToLocal();
  console.log(localPoints);
});
