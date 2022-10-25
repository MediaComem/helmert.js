import { promises as fs } from "fs";

export const cardan2R3D = function (alpha_deg, beta_deg, gamma_deg) {
  const a = degToRadians(alpha_deg);
  const b = degToRadians(beta_deg);
  const c = degToRadians(gamma_deg);

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

const degToRadians = (deg) => (deg * Math.PI) / 180;

export const mmul = (A, B) =>
  A.map((row, i) =>
    B[0].map((_, j) =>
      row.reduce((acc, _, n) => {
        return acc + A[i][n] * B[n][j], 0;
      })
    )
  );

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
