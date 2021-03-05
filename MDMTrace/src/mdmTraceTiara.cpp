#include <cmath>
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

	std::ofstream outputFile;
	std::vector<double> scatteredAngles, beamEnergies, beamPosition;
	double dipoleField=-1, multipoleField=-1;
	//Read from config file
	bool useKinematics = false;
	for(Json::Value::iterator it = config.begin();it!=config.end();it++) {
		if (false) { // DUMMY - AESTHETICS ONLY
		} else if (it.key().asString() == "mdmAngle") {       // MDM ANGLE [degrees]
			mdm->SetMDMAngle(it->asDouble());
			printf("SET: %20s -- %.3f\n","MDM Angle [deg]",mdm->GetMDMAngle());
		} else if (it.key().asString() == "mdmDipoleField") { // MDM FIELD [G]
			// mdm->SetMDMDipoleField(it->asDouble());
			dipoleField = it->asDouble();//.push_back( it->asDouble() );
			printf("SET: %20s -- %.3f\n","MDM Dipole Field [G]",dipoleField);
		} else if(it.key().asString() == "mdmEntranceMultipoleField") { // MULTIPOLE FIELD
			multipoleField = it->asDouble();
			printf("SET: %20s -- %.3f\n","MDM Dipole Field [G]",multipoleField);
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
		} else if(false && it.key().asString() == "beamEnergy") {
			mdm->SetBeamEnergy(it->asDouble());
			printf("SET: %20s -- %.3f\n","Beam Energy [MeV]",mdm->GetBeamEnergy());
    } else if(it.key().asString() == "qValue") {
      mdm->SetQValue(it->asDouble());
      printf("SET: %20s -- %.3f\n","Q Value [MeV]",mdm->GetQValue());
    } else if(it.key().asString() == "excitationEnergy") {
      mdm->SetResidualEnergy(it->asDouble());
      printf("SET: %20s -- %.3f\n","Excitation Energy [MeV]",mdm->GetResidualEnergy());
		}
		else if(it.key().asString() == "scatteredEnergy" ||  // ION ENERGY [MeV]
						it.key().asString() == "scatteredEnergies") {
			if(it->isArray() == false) { // single energy
				mdm->SetScatteredEnergy(it->asDouble());
				printf("SET: %20s -- %.3f\n","Scattered Energy [MeV]",mdm->GetScatteredEnergy());
			}
			else { // multiple energies
				printf("SET: %20s -- ","Scattered Energy [MeV]");
				for(unsigned int i = 0;i<it->size();i++) {
					beamEnergies.push_back((*it)[i].asDouble());
					printf("%.3f ",beamEnergies.back());
				} printf("\n");
			}
		}
#if 0
		else if(it.key().asString() == "scatteredEnergies") { // ION ANGLES [deg]
      for(unsigned int i = 0;i<it->size();i++) {
				beamEnergies.push_back((*it)[i].asDouble());
      }
      printf("SET: %20s -- ","Scattered Angles [deg]");
      for(unsigned int i = 0;i<scatteredAngles.size();i++) {
				printf("%.3f ",scatteredAngles[i]);
      } printf("\n");
		}
#endif
		else if(it.key().asString() == "scatteredAngles") { // ION ANGLES [deg]
      for(unsigned int i = 0;i<it->size();i++) {
				scatteredAngles.push_back((*it)[i].asDouble());
      }
      printf("SET: %20s -- ","Scattered Angles [deg]");
      for(unsigned int i = 0;i<scatteredAngles.size();i++) {
				printf("%.3f ",scatteredAngles[i]);
      } printf("\n");
		}
		else if(it.key().asString() == "beamPosition") { // BEAM POSITION @TARGET [cm]
			if(it->size() == 3){
				for(unsigned int i=0; i< 3; ++i){
					beamPosition.push_back((*it)[i].asDouble());
				}
				printf("SET: %20s -- ", "Beam Position [cm]");
				for(const auto& p : beamPosition){
					printf("%.3f ",p);
				} printf("\n");
			}
		}
		else if(it.key().asString() == "useKinematics") {
      useKinematics = it->asBool();
      if(useKinematics) {
				printf("Calling MDMTrace with kinematics...\n");
      } else {
				printf("Calling MDMTrace without kinematics...\n");
      }
    }
		else if(it.key().asString() == "outputFile") {
			printf("Saving output to file: \"%s\"\n", it->asString().c_str());
			outputFile.open(it->asString().c_str(), std::ofstream::out);
			outputFile << "anglein:energy:angleout:x1:x2:x3:x4/D" << std::endl;
		}
		else if(it.key().asString() != "rayinFile") {
			fprintf(stderr, "WARNING: skipping unrecognized JSON field: %s\n",
							it.key().asString().c_str());
		}
	}
	// Set field values
	if(dipoleField < 0){
		std::cerr << "ERROR: Dipole field not set!\n";
		exit(1);
	}
	if(multipoleField < 0){
		std::cout << "Using standard scaling for entrance multipole field...\n";
		mdm->SetMDMDipoleField(dipoleField);
	} else {
		std::cout << "Using manual value for entrance multipole field...\n";
		mdm->SetMDMDipoleMultipoleField(dipoleField, multipoleField);
	}
	if(beamPosition.empty()){
		mdm->SetBeamPosition(0,0,0);
	} else {
		mdm->SetBeamPosition(
			beamPosition.at(0),beamPosition.at(1),beamPosition.at(2));
	}
	
	auto get_and_print =
		[&](double angle){
			double x1,x2,x3,x4,a1;
			mdm->GetOxfordWirePositions(a1,x1,x2,x3,x4);
			double MM_R1_distance = 38.325 - 2.5;// W1 -> MM Row 1
			double mm1 = x1 + tan(1e-3*a1)*MM_R1_distance;
			//
			double energy = (useKinematics) ? mdm->GetEnergyAfterKinematics() :
				mdm->GetScatteredEnergy();
			printf("Initial Angle/deg: %4.2f  Energy/MeV: %4.2f  W1Angle/mrad (deg):"
						 " %+6.2f (%+3.2f)  Wire Pos/cm: 1: %+5.4f 2: %+5.4f 3: %+5.4f 4: %+5.4f  MM1/cm: %+5.4f\n",
						 angle,energy,a1,a1/17.453293,x1,x2,x3,x4,mm1);
			if(outputFile.is_open()){
				char buf[4096];
				sprintf(buf, "%8.2f\t %8.2f\t %8.2f\t %5.4f\t %5.4f\t %5.4f\t %5.4f",
								angle,energy,a1,x1,x2,x3,x4);
				outputFile << buf << std::endl;
			}
		};

	if(beamEnergies.empty()) {
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
			get_and_print(angle);
		}
	} else {
		for (const auto& energy : beamEnergies) {
			mdm->SetScatteredEnergy(energy);
			for (const auto& angle : scatteredAngles) {
				mdm->SetScatteredAngle(angle);
				mdm->SendRay();
				get_and_print(angle);
			}
		}
	}
	return 0;
}


/////////////////////////////
///// SHUYA'S OLD CODE //////
/////////////////////////////

#if 0		
  for(int i=0;i<11;i++) {
    mdm->SetScatteredAngle(angles[i]);  //ALWAYS
    //mdm->SendRay();                     //WITHOUT KINEMATICS ONLY
    mdm->SendRayWithKinematics();     //KINEMATICS ONLY *
    double x1,x2,x3,x4,a1;
//by Shuya 160712
    mdm->GetOxfordWirePositions(a1,x1,x2,x3,x4);
    //mdm->GetOxfordWirePositions(x1,x2,x3,x4);
//by Shuya 160712
    printf("Initial Angle: %8.2f\t W1Angle: %8.2f\t 1: %5.4f\t 2: %5.4f\t 3: %5.4f\t 4: %5.4f\n",angles[i],a1,x1,x2,x3,x4);
    //printf("Initial Angle: %8.2f\t 1: %5.4f\t 2: %5.4f\t 3: %5.4f\t 4: %5.4f\n",angles[i],x1,x2,x3,x4);
  }
  return 0;
}

  //Set MDM
//  double Brho =  727.7; // kG cm
  //mdm->SetMDMAngle(              0.00);  //ALWAYS SET
  mdm->SetMDMAngle(              0.00);  //ALWAYS SET
//by Shuya 160713. Inputs are in Gauss converted from Tm by 1000.0 * 1000.0 / 160.0.
  //mdm->SetMDMDipoleField(3409.5);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5562.5);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5523.75);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5830.75);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
                                            // (In mdmTrace, ~LISE brho in TM * 1000.0 * 1000.0/160) 
  //mdm->SetMDMDipoleField(5711.5575);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5601.371875);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5523.75);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5417.1875);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(6274.9375);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(6488.285375);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034

//by Shuya 160823. 23Na(d,p)
  //mdm->SetMDMDipoleField(6111.125);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(6051.625);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(6066.625);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5998.375);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(6015.375);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5932.1875);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//23Na (8+), 230 MeV. Brho*1.034
  //mdm->SetMDMDipoleField(8200.2404);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//23Na (11+), 230 MeV. Brho*1.034
  //mdm->SetMDMDipoleField(5963.595);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//24Na (11+), 231.48 MeV. Brho*1.034
  //mdm->SetMDMDipoleField(6111.2502);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(6060.274);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//23Na (11+) beam centered on Micromegas cl3 and 4.
  //mdm->SetMDMDipoleField(5933);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//23Na (11+) beam centered on Micromegas cl4 and 5.
  //mdm->SetMDMDipoleField(5872.086);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//23Na (11+) beam centered on Micromegas cl5 and 6.
  //mdm->SetMDMDipoleField(5809.012);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//23Na (11+) beam centered on Micromegas cl6 and 7.
  //mdm->SetMDMDipoleField(5744.904);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//23Na (11+) beam centered on Micromegas cl2 and 3.
  //mdm->SetMDMDipoleField(5994.615);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//23Na (11+) beam centered on Micromegas cl1 and 2.
  //mdm->SetMDMDipoleField(6054.07);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034

  //mdm->SetMDMDipoleField(5963.625);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034

//by Shuya 160824. 19F(d,p)
  //mdm->SetMDMDipoleField(6198.8125);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5995.);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034

//by Shuya 161003. 19F(d,d) (no scattered beam)
  //mdm->SetMDMDipoleField(6140.926);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(6093.362);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(6033.39);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5971.867);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5906.725);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5841.583);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5781.611);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034

//19F (9+) beam centered on Micromegas cl4.
  //mdm->SetMDMDipoleField(5973.4);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//19F (9+) beam centered on Micromegas cl5.
  //mdm->SetMDMDipoleField(5907.8);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//19F (9+) beam centered on Micromegas cl6.
  //mdm->SetMDMDipoleField(5843.7);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//19F (9+) beam centered on Micromegas cl7.
  //mdm->SetMDMDipoleField(5783.7);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//19F (9+) beam centered on Micromegas cl3.
  //mdm->SetMDMDipoleField(6035.4);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//19F (9+) beam centered on Micromegas cl2.
  //mdm->SetMDMDipoleField(6095.4);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//19F (9+) beam centered on Micromegas cl1.
  //mdm->SetMDMDipoleField(6142.9);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//20F (9+) used for experiment.
  //mdm->SetMDMDipoleField(6198.83);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034


//25Mg (12+) used for experiment.
  //mdm->SetMDMDipoleField(5941.67);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField(5746.3*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (10+) used for experiment.
  //mdm->SetMDMDipoleField(7130.1538);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (9+) used for experiment.
  //mdm->SetMDMDipoleField(7661.9*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (8+) used for experiment.
  //mdm->SetMDMDipoleField(8619.8*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//26Mg (12+) used for experiment.
  //mdm->SetMDMDipoleField(5906.1*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl4-5.
  //mdm->SetMDMDipoleField(5647*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl5-6.
  //mdm->SetMDMDipoleField(5586*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl6-7.
  //mdm->SetMDMDipoleField(5524*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl3-4.
  //mdm->SetMDMDipoleField(5707*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl2-3.
  //mdm->SetMDMDipoleField(5765*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl1-2.
  //mdm->SetMDMDipoleField(5821*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl1-2 adjusted to fit in the data.
  //mdm->SetMDMDipoleField(5906*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl2-3 adjusted to fit in the data.
  //mdm->SetMDMDipoleField(5846*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl3-4 adjusted to fit in the data.
  //mdm->SetMDMDipoleField(5786*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl4-5 adjusted to fit in the data.
  //mdm->SetMDMDipoleField(5723*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl5-6 adjusted to fit in the data.
  //mdm->SetMDMDipoleField(5659*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//25Mg (12+) beam centered on border of Micromegas cl6-7 adjusted to fit in the data.
  //mdm->SetMDMDipoleField(5595*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034


//22Ne(6Li,d). 26Mg(12+) 208.74MeV used for experiment.
  //mdm->SetMDMDipoleField(5351.6*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034

  //mdm->SetMDMDipoleField(3772.7*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//22Ne(6Li,d). 26Mg(12+) 8MeV/u.
  //mdm->SetMDMDipoleField(4783.0*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//22Ne(6Li,d). 26Mg(12+) 7MeV/u.
  //mdm->SetMDMDipoleField(4470.9*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//22Ne(6Li,d).Beam 7MeV/u on the border of micromegas colum 3 an 4.
  //mdm->SetMDMDipoleField(5108.0*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//22Ne(6Li,d).Beam 7MeV/u on the border of micromegas colum 2 an 3.
  //mdm->SetMDMDipoleField(5164.0*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//22Ne(6Li,d).Beam 7MeV/u on the border of micromegas colum 1 an 2.
  //mdm->SetMDMDipoleField(5220.0*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//22Ne(6Li,d).Beam 7MeV/u on the border of micromegas colum 4 an 5.
  //mdm->SetMDMDipoleField(5052.0*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//22Ne(6Li,d).Beam 7MeV/u on the border of micromegas colum 5 an 6.
  //mdm->SetMDMDipoleField(4994.0*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//22Ne(6Li,d).Beam 7MeV/u on the border of micromegas colum 6 an 7.
  //mdm->SetMDMDipoleField(4935.0*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//22Ne(6Li,t).Beam 7MeV/u, 25Mg GS.
  //mdm->SetMDMDipoleField(4465.8*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField((4347.0+30)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField((4589.5-30)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField((4589.5)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034


//59Fe(d,p).Beam 8MeV/u, 60Fe GS.
  //mdm->SetMDMDipoleField((5654.4)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034

//63Cu(d,p).Beam 10MeV/u after 500 ug/cm2 CD2, 63Cu(29+) GS.
  //mdm->SetMDMDipoleField((5928.7355)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//63Cu(d,p).Beam 10MeV/u after 500 ug/cm2 CD2, 63Cu(28+) GS.
  //mdm->SetMDMDipoleField((6140.47389)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//63Cu(d,p).Beam 10MeV/u after 500 ug/cm2 CD2, 63Cu(27+) GS.
  //mdm->SetMDMDipoleField((6367.92795)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//63Cu(d,p).Beam 10MeV/u after 500 ug/cm2 CD2, 63Cu(26+) GS.
  //mdm->SetMDMDipoleField((6612.91103)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//63Cu(d,p).Beam 10MeV/u after 500 ug/cm2 CD2, 63Cu(25+) GS.
  //mdm->SetMDMDipoleField((6877.41779)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  //mdm->SetMDMDipoleField((14477.8)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
  mdm->SetMDMDipoleField((6367.9)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034

//63Cu(d,p).Beam 10MeV/u after 500 ug/cm2 CD2, 64Cu(28+) GS.
  //mdm->SetMDMDipoleField((6198.92408)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034
//63Cu(d,p).Beam 10MeV/u after 500 ug/cm2 CD2, 64Cu(27+) GS.
  //mdm->SetMDMDipoleField((6428.55416)*1.034);  //ALWAYS SET // Brho*1000.0/160. //dipole used*1.034

  //Set Particles                  
//for D(23Na,24Na)p
/*
  mdm->SetTargetMass(    2.00);  //KINEMATICS ONLY *
  mdm->SetProjectileMass( 23.);  //KINEMATICS ONLY *
  mdm->SetScatteredMass(   24.);    //ALWAYS
  //mdm->SetScatteredMass(   23.);    //ALWAYS
  mdm->SetScatteredCharge( 10.);    //ALWAYS
  //Set Energetics
  mdm->SetBeamEnergy(   230.00);  //KINEMATICS ONLY *
  //mdm->SetQValue(        -6.47);  //KINEMATICS ONLY * //Ground state
  mdm->SetQValue(        4.73);  //KINEMATICS ONLY * //Ground state
  //mdm->SetQValue(        0.00);  //KINEMATICS ONLY * //Ground state
  mdm->SetResidualEnergy(0.00);  //KINEMATICS ONLY *
  mdm->SetScatteredEnergy(231.48);  //WITHOUT KINEMATICS ONLY
*/

  //Set Particles                  
//for D(63Cu,64Cu)p
  mdm->SetTargetMass(    2.);  //KINEMATICS ONLY *
  mdm->SetProjectileMass( 63.);  //KINEMATICS ONLY *
  mdm->SetScatteredMass(   63.);    //ALWAYS
  mdm->SetScatteredCharge( 27.);    //ALWAYS
  //mdm->SetBeamEnergy(   619.17);  //KINEMATICS ONLY *
  mdm->SetBeamEnergy(   617.17);  //KINEMATICS ONLY *
  //mdm->SetBeamEnergy(   630.);  //KINEMATICS ONLY *
  //mdm->SetQValue(        5.69);  //KINEMATICS ONLY * //Ground state
  mdm->SetQValue(        0.0);  //KINEMATICS ONLY * //Ground state
  mdm->SetResidualEnergy(0.0);  //KINEMATICS ONLY *

  //Set Energetics
  //mdm->SetBeamEnergy(   472.);  //KINEMATICS ONLY *
  //mdm->SetQValue(        4.38);  //KINEMATICS ONLY * //Ground state
  //mdm->SetQValue(        6.87);  //KINEMATICS ONLY * //Ground state
  ///mdm->SetQValue(        9.14);  //KINEMATICS ONLY * //Ground state
  //mdm->SetQValue(        8.87);  //KINEMATICS ONLY * //Ground state
  //mdm->SetQValue(        0.0);  //KINEMATICS ONLY * //Ground state
  //mdm->SetQValue(        8.15);  //KINEMATICS ONLY * //Ground state
  //mdm->SetQValue(        4.30);  //KINEMATICS ONLY * //Ground state
  //mdm->SetQValue(        9.14);  //KINEMATICS ONLY * //Ground state
  //mdm->SetResidualEnergy(11.10);  //KINEMATICS ONLY *
  //mdm->SetResidualEnergy(0.0);  //KINEMATICS ONLY *
	//20F GS
  //mdm->SetScatteredEnergy(191.26);  //WITHOUT KINEMATICS ONLY
	//20F Ex10MeV
  //mdm->SetScatteredEnergy(184.03);  //WITHOUT KINEMATICS ONLY
	//18O GS
  //mdm->SetScatteredEnergy(185.81);  //WITHOUT KINEMATICS ONLY
	//18O Ex10MeV
  //mdm->SetScatteredEnergy(169.85);  //WITHOUT KINEMATICS ONLY
	//17O GS
  //mdm->SetScatteredEnergy(199.13);  //WITHOUT KINEMATICS ONLY
	//17O Ex10MeV
  //mdm->SetScatteredEnergy(186.45);  //WITHOUT KINEMATICS ONLY
	//19F GS
  //mdm->SetScatteredEnergy(190.00);  //WITHOUT KINEMATICS ONLY
	//19F Ex10MeV
  //mdm->SetScatteredEnergy(178.20);  //WITHOUT KINEMATICS ONLY
	//26Mg GS
  //mdm->SetScatteredEnergy(110);  //WITHOUT KINEMATICS ONLY
  //mdm->SetScatteredEnergy(101.96);  //WITHOUT KINEMATICS ONLY
  //mdm->SetScatteredEnergy(250.0);  //WITHOUT KINEMATICS ONLY

//25Mg (from 22Ne 8MeV/n->26Mg-neutron->25Mg, 1deg)
  //mdm->SetScatteredEnergy(160.473);  //WITHOUT KINEMATICS ONLY
//25Mg (from 22Ne 8MeV/n->26Mg-neutron->25Mg, 0deg)
  //mdm->SetScatteredEnergy(160.81);  //WITHOUT KINEMATICS ONLY
//25Mg (from 22Ne 8MeV/n->26Mg-neutron->25Mg, 2deg)
  //mdm->SetScatteredEnergy(159.4);  //WITHOUT KINEMATICS ONLY
//25Mg (from 22Ne 7MeV/n->26Mg-neutron->25Mg, 1deg)
  //mdm->SetScatteredEnergy(140.27);  //WITHOUT KINEMATICS ONLY
//25Mg (from 22Ne 8MeV/n->26Mg-neutron->25Mg, 0deg)
  //mdm->SetScatteredEnergy(140.58);  //WITHOUT KINEMATICS ONLY
//25Mg (from 22Ne 8MeV/n->26Mg-neutron->25Mg, 2deg)
  //mdm->SetScatteredEnergy(139.35);  //WITHOUT KINEMATICS ONLY

  //mdm->SetScatteredEnergy(243.61);  //WITHOUT KINEMATICS ONLY
  //mdm->SetScatteredEnergy(249.04);  //WITHOUT KINEMATICS ONLY
  //mdm->SetScatteredEnergy(254.00);  //WITHOUT KINEMATICS ONLY
  //mdm->SetScatteredEnergy(153.6);  //WITHOUT KINEMATICS ONLY
  //mdm->SetScatteredEnergy(140);  //WITHOUT KINEMATICS ONLY
  //mdm->SetScatteredEnergy(137.5);  //WITHOUT KINEMATICS ONLY
  //mdm->SetScatteredEnergy(142.5);  //WITHOUT KINEMATICS ONLY



/*
//for 6Li(22Ne,26Mg)d
  mdm->SetTargetMass(    6.00);  //KINEMATICS ONLY *
  mdm->SetProjectileMass( 22.);  //KINEMATICS ONLY *
  mdm->SetScatteredMass(   26.);    //ALWAYS
  mdm->SetScatteredCharge( 12.);    //ALWAYS
  //Set Energetics
  mdm->SetQValue(        9.14);  //KINEMATICS ONLY * //Ground state
  mdm->SetResidualEnergy(10.2);  //KINEMATICS ONLY *
  //if 10MeV/n
  //mdm->SetBeamEnergy(   220.00);  //KINEMATICS ONLY *
  ///mdm->SetScatteredEnergy(208.0);  //WITHOUT KINEMATICS ONLY
  //if 8MeV/n
  mdm->SetBeamEnergy(   186.00);  //KINEMATICS ONLY *
  mdm->SetScatteredEnergy(176.3);  //WITHOUT KINEMATICS ONLY
  //if 5MeV/n
  //mdm->SetBeamEnergy(   110.00);  //KINEMATICS ONLY *
  //mdm->SetScatteredEnergy(104.);  //WITHOUT KINEMATICS ONLY
*/

/*for 6Li(22Ne,22Ne)6Li Elastic
  mdm->SetTargetMass(    6.00);  //KINEMATICS ONLY *
  mdm->SetProjectileMass( 22.);  //KINEMATICS ONLY *
  mdm->SetScatteredMass(   22.);    //ALWAYS
  mdm->SetScatteredCharge( 10.);    //ALWAYS
  //Set Energetics
  mdm->SetBeamEnergy(   220.00);  //KINEMATICS ONLY *
  mdm->SetQValue(        0.0);  //KINEMATICS ONLY * //Ground state
  mdm->SetResidualEnergy(0.0);  //KINEMATICS ONLY *
  mdm->SetScatteredEnergy(208.);  //WITHOUT KINEMATICS ONLY
*/

  //Set Angles
  double angle = 0.00;
  //double angles[7] =  {1.95,3.42,4.185,4.95,5.715,6.48,7.95};
  //double angles[7] = {angle-(2.0*0.765 + 1.47),angle-(2.0*0.765),angle-0.765,angle,angle+0.765,angle+(2.0*0.765),angle + (2.0*0.765 + 1.47)};
  //double angles[11] = {-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5};
//Only for heavy recoils such as Cu which are not scattered above 1.5 deg in Lab sys.
  double angles[11] = {-1.45,-1.2,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.2,1.45};

  printf("Oxford Wire Positions\n");
//by Shuya 160713.
  //for(int i=0;i<7;i++) {
  for(int i=0;i<11;i++) {
    mdm->SetScatteredAngle(angles[i]);  //ALWAYS
    //mdm->SendRay();                     //WITHOUT KINEMATICS ONLY
    mdm->SendRayWithKinematics();     //KINEMATICS ONLY *
    double x1,x2,x3,x4,a1;
//by Shuya 160712
    mdm->GetOxfordWirePositions(a1,x1,x2,x3,x4);
    //mdm->GetOxfordWirePositions(x1,x2,x3,x4);
//by Shuya 160712
    printf("Initial Angle: %8.2f\t W1Angle: %8.2f\t 1: %5.4f\t 2: %5.4f\t 3: %5.4f\t 4: %5.4f\n",angles[i],a1,x1,x2,x3,x4);
    //printf("Initial Angle: %8.2f\t 1: %5.4f\t 2: %5.4f\t 3: %5.4f\t 4: %5.4f\n",angles[i],x1,x2,x3,x4);
  }
  return 0;
}

// positions in cm
// angles in mrad
#endif
