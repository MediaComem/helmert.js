import { describe, expect, test } from "@jest/globals";
import { Cardan2R3D } from "../lib/helmert.js";

describe("A rotation matrix from Cardan Angles", () => {
  const returnedMatrix = Cardan2R3D(1, 1, 1);

  test("is 3x3", () => {
    expect(returnedMatrix).toHaveLength(3);
    returnedMatrix.forEach((subarray) => {
      expect(subarray).toHaveLength(3);
    });
  });

  test("has expected values", () => {
    const expectedValues = [
      0.99969541, 0.01775429, -0.0171425, -0.01744975, 0.9996901, 0.01775429,
      0.01745241, -0.01744975, 0.99969541
    ].flat();
    const flatMatrix = returnedMatrix.flat();

    flatMatrix.forEach((num, index) => {
      expect(num).toBeCloseTo(expectedValues[index], 8);
    });
  });
});
