import random
import helmertTransformation

def generateGlobalPoint():
    x = random.uniform(4346300, 4346450)
    y = random.uniform(507400, 507550)
    z = random.uniform(4625200, 4625400)
    return [x, y, z]


def generateLocalPoint():
    x = random.uniform(-100, 0)
    y = random.uniform(-100, 100)
    z = random.uniform(0, 150)
    return [x, y, z]


def generateFile(fn, num_points, path):
    file_out = open(path, "w")
    test_data = []
    for i in range(num_points):
        data = fn()
        key = "P"+str(i)
        file_out.write('{}\t{:0.6f}\t{:0.6f}\t{:0.6f}{newline}'.format(
            key, data[0], data[1], data[2], newline="\n" if i != num_points - 1 else ""))
    return test_data


def generateTestData(num_points):
    generateFile(generateGlobalPoint, num_points,
                 "./data/generated_data/testGlobalData.xyz")
    generateFile(generateLocalPoint, num_points,
                 "./data/generated_data/testLocalData.xyz")


def testHelmert(num_points):
    print("GENERATING DATA FOR TESTING HELMERT ALGORITHM")
    generateTestData(num_points)
    helm = helmertTransformation.HelmertTransformation3D('1')
    helm.importPointsGlobal('./data/generated_data/testGlobalData.xyz')
    helm.importPointsLocal('./data/generated_data/testLocalData.xyz')
    helm.estimateHelmert3DSVD_minimum()
    points_local = helmertTransformation.global2local(helm, helm.points_global)
    print("WRITING HELMERT RESULTS IN PY/DATA/GENERATED_DATA/RESULT.XYZ")
    file_out = open("./data/generated_data/result.xyz", "w")
    for index, key in enumerate(points_local):
        file_out.write('{}\t{:0.6f}\t{:0.6f}\t{:0.6f}{newline}'.format(
            key, points_local[key][0], points_local[key][1], points_local[key][2], newline="\n" if index != num_points - 1 else ""))
    

testHelmert(10)
