###################################################################################################
EXE=DRA
test: clear mesh run view
mesh: bmesh
###############
clean:
	(rm -rf *.log [1-9]*);
clear:
	(./Allclean)
###############
bmesh:
#	(cp -rf ./system/controlDict.mesh ./system/controlDict);
	(blockMesh);
regionMesh:
#air region
	(rm -rf constant/polyMesh/sets constant/polyMesh/*Zones);
	(setSet -batch ./config/make.setSetAir -noVTK);
	(rm -rf constant/polyMesh/sets/*_old);
	(subsetMesh air);
	(cp -rf 1/* ./constant/air/);
	(rm -rf 1);

#fuel region
	(rm -rf constant/polyMesh/sets constant/polyMesh/*Zones);
	(setSet -batch ./config/make.setSetFuel -noVTK);
	(rm -rf constant/polyMesh/sets/*_old);
	(subsetMesh fuel);
	(cp -rf 1/* ./constant/fuel/);
	(rm -rf 1);
#interconnect region
	(rm -rf constant/polyMesh/sets constant/polyMesh/*Zones);
	(setSet -batch ./config/make.setSetInterconnect -noVTK);
	(rm -rf constant/polyMesh/sets/*_old);
	(subsetMesh interconnect);
	(cp -rf 1/* ./constant/interconnect/);
	(rm -rf 1);
#electrolyte region
	(rm -rf constant/polyMesh/sets constant/polyMesh/*Zones);
	(setSet -batch ./config/make.setSetElectrolyte -noVTK);
	(rm -rf constant/polyMesh/sets/*_old);
	(subsetMesh electrolyte);
	(cp -rf 1/* ./constant/electrolyte/);
	(rm -rf 1);

#redefine boundary conditions air region
	(rm -rf ./constant/air/polyMesh/sets);
	(setSet -batch ./config/make.faceSetAir -region air -noVTK);
	(rm -rf ./constant/air/polyMesh/sets/*_old);
	(setSet -batch ./config/make.faceAir -region air -noVTK);
	(rm -rf ./constant/air/polyMesh/sets/*_old);

	(setsToZones -noFlipMap -region air);
	(createPatch -overwrite -region air);

#redifine boundary conditions fuel region
	(rm -rf ./constant/fuel/polyMesh/sets);
	(setSet -batch ./config/make.faceSetFuel -region fuel -noVTK);
	(rm -rf ./constant/fuel/polyMesh/sets/*_old);
	(setSet -batch ./config/make.faceFuel -region fuel -noVTK);
	(rm -rf ./constant/fuel/polyMesh/sets/*_old);

	(setsToZones -noFlipMap -region fuel);
	(createPatch -overwrite -region fuel);

#zones for mapping air
	(rm -rf constant/air/polyMesh/sets constant/air/polyMesh/*Zones);
	(setSet -batch ./config/make.setAir -region air -noVTK);

	(setSet -batch ./config/make.faceSetTest -region air -noVTK);
	(setSet -batch ./config/make.faceTest -region air -noVTK);
	(rm -rf constant/air/polyMesh/sets/*_old);

	(setsToZones -noFlipMap -region air);

#zones for mapping fuel
	(rm -rf constant/fuel/polyMesh/sets constant/fuel/polyMesh/*Zones);
	(setSet -batch ./config/make.setFuel -region fuel -noVTK);

	(setSet -batch ./config/make.faceSetTest -region fuel -noVTK);
	(setSet -batch ./config/make.faceTest -region fuel -noVTK);
	(rm -rf constant/fuel/polyMesh/sets/*_old);

	(setsToZones -noFlipMap -region fuel);

#zones for mapping interconnect
	(rm -rf constant/interconnect/polyMesh/sets constant/interconnect/polyMesh/*Zones);
	(setSet -batch ./config/make.setInterconnect -region interconnect -noVTK);
	(rm -rf constant/interconnect/polyMesh/sets/*_old);
	(setsToZones -noFlipMap -region interconnect);
#zones for mapping electrolyte
	(rm -rf constant/electrolyte/polyMesh/sets constant/electrolyte/polyMesh/*Zones);
	(setSet -batch ./config/make.setElectrolyte -region electrolyte -noVTK);
	(rm -rf constant/electrolyte/polyMesh/sets/*_old);
	(setsToZones -noFlipMap -region electrolyte);


patchBC:

control:
#	(cp -rf ./system/controlDict.run ./system/controlDict);
#############
run:
	(make clean);
	(electricFoamAxial | tee srun.log);
post:
	(sample -region air);
view:
	(foamToVTK -latestTime -surfaceFields);
	(foamToVTK -latestTime -cellSet nonIrradiated);
#	(foamToVTK -latestTime -region air);
#	(foamToVTK -latestTime -region fuel);
#	(foamToVTK -latestTime -region interconnect);
#	(foamToVTK -latestTime -region electrolyte);

viewAll:
	(foamToVTK);
	(foamToVTK -region air);
	(foamToVTK -region fuel);
	(foamToVTK -region electrolyte);
	(foamToVTK -region interconnect0);
	(foamToVTK -region interconnect1);
###################################################################################################
