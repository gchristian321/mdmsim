#include <stdio.h>
#include <iostream>
#include <fstream>
#include "json/json.h"
#include "MDMTrace.h"
#include "Rayin.h"


// Example program to propagate a defined ray through
// the MDM.
//
// Reads inputs from a json config file. See the example
// config-mdmExample.json 
//
int main(int argc, char* argv[]) {

	// Check for correct arguments
  if(argc<2) {
    std::cout << "Usage: mdmTrace config-file" << std::endl;
    return 0;
  }

  Json::Value config;
  std::string configFileName = argv[1];
  std::ifstream configStream(configFileName.c_str());
  configStream >> config;
  configStream.close();

	// Create link to rayin.dat file
	// NOTE:: Must come before the first call to
	// MDMTrace::Instance()
	Json::Value::iterator itRay = config.begin();
	while(itRay != config.end()) {
		if (itRay.key().asString() == "rayinFile"){ // RAYTRACE FILE
			break;
		}
		++itRay;
	}
	if(itRay == config.end()) {
		std::cerr << "ERROR: \"rayinFile\" not set in \"" << argv[1] << "\"\n";
		exit(1);
	}
	Rayin rayin_(itRay->asString());

	
  //Get Global MDM Instance
  MDMTrace* mdm = MDMTrace::Instance();

	std::vector<double> scatteredAngles;
	//Read from config file
	bool useKinematics = false;
	for(Json::Value::iterator it = config.begin();it!=config.end();it++) {
		if (false) { // DUMMY - AESTHETICS ONLY
		}
		else if (it.key().asString() == "mdmAngle") {       // MDM ANGLE [degrees]
			mdm->SetMDMAngle(it->asDouble());
			printf("SET: %20s -- %.3f\n","MDM Angle [deg]",mdm->GetMDMAngle());
		}
		else if (it.key().asString() == "mdmDipoleField") { // MDM FIELD [G]
			mdm->SetMDMDipoleField(it->asDouble());
			printf("SET: %20s -- %.3f\n","MDM Dipole Field [G]",mdm->GetMDMDipoleField());
		} else if(it.key().asString() == "targetMass") {
      mdm->SetTargetMass(it->asDouble());
      printf("SET: %20s -- %.3f\n","Target Mass [amu]",mdm->GetTargetMass());
    } else if(it.key().asString() == "projectileMass") {
      mdm->SetProjectileMass(it->asDouble());
      printf("SET: %20s -- %.3f\n","Projectile Mass [amu]",mdm->GetProjectileMass());
		} else if (it.key().asString() == "scatteredMass") {  // ION MASS [AMU]
      mdm->SetScatteredMass(it->asDouble());
      printf("SET: %20s -- %.3f\n","Scattered Mass [AMU]",mdm->GetScatteredMass());
		}
		else if(it.key().asString() == "scatteredCharge") { // ION CHARGE STATE [e+ charge]
			mdm->SetScatteredCharge(it->asDouble());
      printf("SET: %20s -- %.3f\n","Scattered Charge [e+]",mdm->GetScatteredCharge());
		} else if(it.key().asString() == "beamEnergy") {
      mdm->SetBeamEnergy(it->asDouble());
      printf("SET: %20s -- %.3f\n","Beam Energy [MeV]",mdm->GetBeamEnergy());
    } else if(it.key().asString() == "qValue") {
      mdm->SetQValue(it->asDouble());
      printf("SET: %20s -- %.3f\n","Q Value [MeV]",mdm->GetQValue());
    } else if(it.key().asString() == "excitationEnergy") {
      mdm->SetResidualEnergy(it->asDouble());
      printf("SET: %20s -- %.3f\n","Excitation Energy [MeV]",mdm->GetResidualEnergy());
		}	else if(it.key().asString() == "scatteredEnergy") { // ION ENERGY [MeV]
      mdm->SetScatteredEnergy(it->asDouble());
      printf("SET: %20s -- %.3f\n","Scattered Energy [MeV]",mdm->GetScatteredEnergy());
    }
		else if(it.key().asString() == "scatteredAngles") { // ION ANGLES [deg]
      for(unsigned int i = 0;i<it->size();i++) {
				scatteredAngles.push_back((*it)[i].asDouble());
      }
      printf("SET: %20s -- ","Scattered Angles [deg]");
      for(unsigned int i = 0;i<scatteredAngles.size();i++) {
				printf("%.3f ",scatteredAngles[i]);
      } printf("\n");
		} else if(it.key().asString() == "useKinematics") {
      useKinematics = it->asBool();
      if(useKinematics) {
				printf("Calling MDMTrace with kinematics...\n");
      } else {
				printf("Calling MDMTrace without kinematics...\n");
      }
    }
	}
	for (const auto& angle : scatteredAngles) {
		mdm->SetScatteredAngle(angle);
		
		if(useKinematics)
		{
			mdm->SendRayWithKinematics();
		}
		else
		{
			mdm->SendRay();
		}
		
    double x1,x2,x3,x4,a1;
    mdm->GetOxfordWirePositions(a1,x1,x2,x3,x4);
		//
		double energy = (useKinematics) ? mdm->GetEnergyAfterKinematics() : mdm->GetScatteredEnergy();
    printf("Initial Angle: %8.2f\t Energy: %8.2f\t W1Angle: %8.2f\t 1: %5.4f\t 2: %5.4f\t 3: %5.4f\t 4: %5.4f\n",
					 angle,energy,a1,x1,x2,x3,x4);
	}

	return 0;
}
