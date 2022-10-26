import helmertTransformation 


helm = helmertTransformation.HelmertTransformation3D('1')





helm.importPointsGlobal('trajWGS84.xyz')
helm.importPointsLocal('trajLocal.xyz')

# helm.plot2DMap(typePoints= 'local')
# helm.plot2DMap(typePoints= 'global')

# # nb_iterations_ransac = 100
# # helm.estimateHelmert3DSVD(nb_iterations_ransac)
# # helm.printParameters()
# # helm.printResiduals()

helm.estimateHelmert3DSVD_minimum()
# helm.printParameters()
# helm.printResiduals()

# helm.exportAllLocalPointsINGlobalCoord('ptsLocalTransfInWGS84.xyz')

# points_local = helmertTransformation.global2local(helm, helm.points_global)
# points_global = helmertTransformation.local2global(helm, helm.points_local)

# print(points_local)