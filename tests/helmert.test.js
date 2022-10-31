import { describe, expect, test } from "@jest/globals";
import {
  cardan2R3D,
  readXYZFile,
  parseXYZData,
  R3D2Cardan
} from "../src/utils.js";
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
      expect(num).toBeCloseTo(expectedValues[index], 6);
    });
  });
});

describe("Cardan Angles from a rotation Matrix", () => {
  const rotationMatrix = [
    [0.98137454, 0.13181889, -0.13974186],
    [-0.0919062, 0.96094634, 0.26102756],
    [0.16869279, -0.24332266, 0.95516325]
  ];

  const returnedValues = R3D2Cardan(rotationMatrix);

  test("has expected values", () => {
    const expectedValues = [
      14.291809812321672, 9.711823613140364, 5.350172730672064
    ];

    returnedValues.forEach((num, index) => {
      expect(num).toBeCloseTo(expectedValues[index], 6);
    });
  });
});

describe("Helmert algorithms produce same results in JS and Python", () => {
  const helm = new Helmert();
  test("has expected values", async () => {
    const expectedValues = Object.values(
      await readXYZFile("./data/generated_data/result.xyz").then((d) =>
        parseXYZData(d)
      )
    ).flat();
    const actualValues = Object.values(
      await helm
        .importFiles(
          "./data/generated_data/testLocalData.xyz",
          "./data/generated_data/testGlobalData.xyz"
        )
        .then(() => helm.estimateHelmertMinimum().globalToLocal())
    ).flat();
    console.log(actualValues, expectedValues)

    actualValues.forEach((num, index) => {
      expect(num).toBeCloseTo(expectedValues[index], 6);
    });
  });
});
