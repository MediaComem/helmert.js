import { describe, expect, test } from "@jest/globals";
import { cardan2R3D, readXYZFile, parseXYZData } from "../src/utils.js";
import { Helmert } from "../src/helmert.js";

describe("A rotation matrix from Cardan Angles", () => {
  const returnedMatrix = cardan2R3D(1, 1, 1);

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

describe("Helmert and global to local algorithms produce same results in JS and Python", () => {
  const helm = new Helmert();
  test("has expected values", async () => {
    const expectedValues = Object.values(
      await readXYZFile("./data/generated_data/result.xyz").then((d) =>
        parseXYZData(d)
      )
    ).flat();
    const actualValues = Object.values(await helm
      .importFiles(
        "./data/generated_data/testLocalData.xyz",
        "./data/generated_data/testGlobalData.xyz"
      )
      .then(() => helm.estimateHelmertMinimum().globalToLocal())).flat();
      actualValues.forEach((num, index) => {
        expect(num).toBeCloseTo(expectedValues[index], 8)
      })
  });
});
