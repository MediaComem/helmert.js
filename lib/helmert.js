export const Cardan2R3D = function (
  alpha_deg,
  beta_deg,
  gamma_deg
) {
  const a = degToRadians(alpha_deg);
  const b = degToRadians(beta_deg);
  const c  = degToRadians(gamma_deg);

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

const mmul = (A, B) => {
  const result = new Array(A.length)
    .fill(0)
    .map((row) => new Array(B[0].length).fill(0));
  return result.map((row, i) => {
    return row.map((val, j) => {
      return A[i].reduce((sum, elm, k) => sum + elm * B[k][j], 0);
    });
  });
};
