import { eye, readXYZFile, parseXYZData, zeros, mmul } from "./utils.js";
import { SVD } from "svd-js";

export class Helmert {
  constructor() {
    this.successful = false;
    this.iter = "99 (99)";
    this.JD_UTC = NaN;
    this.points_global = {};
    this.points_local = {};
    this.points_common = {};
    this.helmert3DParam = {
      tX: [0.0, "unknown"],
      tY: [0.0, "unknown"],
      tZ: [0.0, "unknown"],
      rX: [0.0, "unknown"],
      rY: [0.0, "unknown"],
      rZ: [0.0, "unknown"]
    };

    this.helmert3DParamPrecision = {
      s0: [NaN, "unknown"],
      stX: [NaN, "unknown"],
      stY: [NaN, "unknown"],
      stZ: [NaN, "unknown"],
      srX: [NaN, "unknown"],
      srY: [NaN, "unknown"],
      srZ: [NaN, "unknown"]
    };

    this.sigma_X = 0.001;
    this.sigma_Y = 0.001;
    this.sigma_Z = 0.001;
    this.sigma_x = 0.001;
    this.sigma_y = 0.001;
    this.sigma_z = 0.001;

    this.sigma = 0.001;
  }
  async importFiles(localPath, globalPath) {
    this.points_global = await readXYZFile(globalPath).then((d) =>
      parseXYZData(d)
    );
    this.points_local = await readXYZFile(localPath).then((d) =>
      parseXYZData(d)
    );
  }
  estimateHelmertMinimum() {
    // Estimate Helmert 3D with singular values decomposition method
    this.points_common = {};
    for (let key in this.points_global) {
      if (this.points_local.hasOwnProperty(key)) {
        let X = this.points_global[key][0];
        let Y = this.points_global[key][1];
        let Z = this.points_global[key][2];
        let x = this.points_local[key][0];
        let y = this.points_local[key][1];
        let z = this.points_local[key][2];
        this.points_common[key] = [X, Y, Z, x, y, z, NaN, NaN, NaN, NaN, false];
      }
    }

    if (Object.keys(this.points_common).length < 3) {
      console.warn("NOT ENOUGH POINTS FOR HELMERT3D");
      this.successful = false;
      return;
    }

    // Create weight matrix
    const P = eye(Object.keys(this.points_common).length * 3);

    // Compute Helmert adjustment with 3 best points
    this.computeHelmert(this.points_common, P);
  }

  computeHelmert(points, weights) {
    // Compute local and global centroid
    let global_X = 0;
    let global_Y = 0;
    let global_Z = 0;
    let local_x = 0;
    let local_y = 0;
    let local_z = 0;
    let poids_X = 0;
    let poids_Y = 0;
    let poids_Z = 0;
    let weight_x = 0;
    let weight_y = 0;
    let weight_z = 0;
    let k = 0;

    for (let key in points) {
      global_X += points[key][0] * weights[k][k];
      global_Y += points[key][1] * weights[k + 1][k + 1];
      global_Z += points[key][2] * weights[k + 2][k + 2];

      poids_X += weights[k][k];
      poids_Y += weights[k + 1][k + 1];
      poids_Z += weights[k + 2][k + 2];

      local_x += points[key][3] * weights[k][k];
      local_y += points[key][4] * weights[k + 1][k + 1];
      local_z += points[key][5] * weights[k + 2][k + 2];

      weight_x += weights[k][k];
      weight_y += weights[k + 1][k + 1];
      weight_z += weights[k + 2][k + 2];
    }
    // Global centroid
    const centroid_X = global_X / poids_X;
    const centroid_Y = global_Y / poids_Y;
    const centroid_Z = global_Z / poids_Z;

    // Local centroid
    const centroid_x = local_x / weight_x;
    const centroid_y = local_y / weight_y;
    const centroid_z = local_z / weight_z;

    const global_centroid = [[centroid_X], [centroid_Y], [centroid_Z]];
    const local_centroid = [[centroid_x], [centroid_y], [centroid_z]];

    const W_SVD = eye(Object.keys(points).length);

    // Compute Centered Vectors
    const global_matrix = zeros(Object.keys(points).length, 3);
    const local_matrix = zeros(Object.keys(points).length, 3);

    k = 0;

    for (let key in points) {
      let reduit_X = points[key][0] - centroid_X;
      let reduit_Y = points[key][1] - centroid_Y;
      let reduit_Z = points[key][2] - centroid_Z;

      let reduit_x = points[key][3] - centroid_x;
      let reduit_y = points[key][4] - centroid_y;
      let reduit_z = points[key][5] - centroid_z;

      global_matrix[0][k] = reduit_X;
      global_matrix[1][k] = reduit_Y;
      global_matrix[2][k] = reduit_Z;

      local_matrix[0][k] = reduit_x;
      local_matrix[1][k] = reduit_y;
      local_matrix[2][k] = reduit_z;
      k += 1;
    }

    // Compute covariance matrix
    
  }
}
