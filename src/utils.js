import { promises as fs } from "fs";

export const cardan2R3D = (alpha_deg, beta_deg, gamma_deg) => {
  const a = (alpha_deg * Math.PI) / 180;
  const b = (beta_deg * Math.PI) / 180;
  const c = (gamma_deg * Math.PI) / 180;

  const Rz = [
    [Math.cos(c), Math.sin(c), 0],
    [-Math.sin(c), Math.cos(c), 0],
    [0, 0, 1]
  ];
  const Ry = [
    [Math.cos(b), 0, -Math.sin(b)],
    [0, 1, 0],
    [Math.sin(b), 0, Math.cos(b)]
  ];
  const Rx = [
    [1, 0, 0],
    [0, Math.cos(a), Math.sin(a)],
    [0, -Math.sin(a), Math.cos(a)]
  ];

  const _r = mmul(Rz, Ry);
  const R = mmul(_r, Rx);
  return R;
};

export const R3D2Cardan = (R) => {
  let a = 0;
  let b = 0;
  let c = 0;
  if (R[2][0] === 1 || R[2][0] === -1) {
    c = 0;
    if (R[2][0] === -1) {
      b = Math.PI / 2;
      a = c + Math.atan2(R[0][1], R[0][2]);
    } else {
      b = -Math.PI / 2;
      a = -c + Math.atan2(R[0][1], R[0][2]);
    }
  } else {
    b = -Math.asin(R[2][0]);
    a = Math.atan2(R[2][1] / Math.cos(b), R[2][2] / Math.cos(b));
    c = Math.atan2(R[1][0] / Math.cos(b), R[0][0] / Math.cos(b));
  }
  const alpha_deg = (-a * 180) / Math.PI;
  const beta_deg = (-b * 180) / Math.PI;
  const gamma_deg = (-c * 180) / Math.PI;
  return [alpha_deg, beta_deg, gamma_deg];
};

export const readXYZFile = async (path) => {
  const data = await fs
    .readFile(path, "utf8")
    .catch((err) => console.error("Failed to read XYZ file", err));
  return data;
};

export const parseXYZData = (data) => {
  const fileData = {};
  for (const line of data.split("\n")) {
    const parsedLine = line.split("\t");
    fileData[parsedLine[0]] = [
      +parsedLine[1].trim(),
      +parsedLine[2].trim(),
      +parsedLine[3].trim()
    ];
  }
  return fileData;
};

export const eye = (N, M = N, k = 0) => {
  const t = [];
  for (let i = 0; i < N; i++) {
    const p = [];
    for (let j = 0; j < M; j++) {
      p.push(j - i === k ? 1 : 0);
    }
    t.push(p);
  }
  return t;
};

export const zeros = (w, h, v = 0) =>
  Array.from(new Array(h), (_) => Array(w).fill(v));

export const mmul = (a, b) => {
  if (!Array.isArray(a) || !Array.isArray(b) || !a.length || !b.length) {
    throw new Error("arguments should be in 2-dimensional array format");
  }
  let x = a.length,
    z = a[0].length,
    y = b[0].length;
  if (b.length !== z) {
    throw new Error(
      "number of columns in the first matrix should be the same as the number of rows in the second"
    );
  }
  let productRow = Array.apply(null, new Array(y)).map(
    Number.prototype.valueOf,
    0
  );
  let product = new Array(x);
  for (let p = 0; p < x; p++) {
    product[p] = productRow.slice();
  }
  for (let i = 0; i < x; i++) {
    for (let j = 0; j < y; j++) {
      for (let k = 0; k < z; k++) {
        product[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  return product;
};

export const mtrans = (matrix) =>
  matrix[0].map((_, colIndex) => matrix.map((row) => row[colIndex]));

export const mdet = (matrix) =>
  matrix.length == 1
    ? matrix[0][0]
    : matrix.length == 2
    ? matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    : matrix[0].reduce(
        (r, e, i) =>
          r +
          (-1) ** (i + 2) *
            e *
            mdet(matrix.slice(1).map((c) => c.filter((_, j) => i != j))),
        0
      );
