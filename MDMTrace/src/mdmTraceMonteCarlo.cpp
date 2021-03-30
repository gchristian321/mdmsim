#include <cmath>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <array>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include "json/json.h"
#include "MDMTrace.h"
#include "Rayin.h"


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
		return 1;
	}
	Rayin rayin_(itRay->asString());

	
  //Get Global MDM Instance
  MDMTrace* mdm = MDMTrace::Instance();

	std::string outputFileName = "mdmTraceMonteCarlo.root";
	std::vector<double> scatteredAngles;
	std::vector<double> dipoleField, multipoleField;
	size_t numberOfRays = 1000;
	bool exIsScattered = true;
	//
	// Read from config file
	for(Json::Value::iterator it = config.begin();it!=config.end();it++) {
		if (it.key().asString() == "mdmAngle") {       // MDM ANGLE [degrees]
			mdm->SetMDMAngle(it->asDouble());
			printf("SET: %20s -- %.3f\n","MDM Angle [deg]",mdm->GetMDMAngle());
		}
		else if (it.key().asString() == "mdmDipoleField" ||  // MDM FIELD [G]
						 it.key().asString() == "mdmDipoleHallProbe")
		{
			// scale field value if hall probe entered instead of real field
			double scl = it.key().asString() == "mdmDipoleField" ? 1. : 1.034;
			dipoleField.push_back( scl*it->asDouble() );
			printf("SET: %20s -- %.3f\n","MDM Dipole Field [G]", dipoleField.back());
		}
		else if(it.key().asString() == "mdmEntranceMultipoleField") { // MULTIPOLE FIELD
			multipoleField.push_back( it->asDouble() );
			printf("SET: %20s -- %.3f\n","Multipole Field [G]",multipoleField.at(0));
		}
		else if(it.key().asString() == "targetMass") {
      mdm->SetTargetMass(it->asDouble());
      printf("SET: %20s -- %.3f\n","Target Mass [amu]",mdm->GetTargetMass());
    }
		else if(it.key().asString() == "projectileMass") {
      mdm->SetProjectileMass(it->asDouble());
      printf("SET: %20s -- %.3f\n","Projectile Mass [amu]",mdm->GetProjectileMass());
		}
		else if (it.key().asString() == "scatteredMass") {  // ION MASS [AMU]
      mdm->SetScatteredMass(it->asDouble());
      printf("SET: %20s -- %.3f\n","Scattered Mass [AMU]",mdm->GetScatteredMass());
		}
		else if(it.key().asString() == "scatteredCharge") { // ION CHARGE STATE [e+ charge]
			mdm->SetScatteredCharge(it->asDouble());
      printf("SET: %20s -- %.3f\n","Scattered Charge [e+]",mdm->GetScatteredCharge());
		}
		else if(it.key().asString() == "beamEnergy") {
			mdm->SetBeamEnergy(it->asDouble());
			printf("SET: %20s -- %.3f\n","Beam Energy [MeV]",mdm->GetBeamEnergy());
    }
		else if(it.key().asString() == "qValue") {
      mdm->SetQValue(it->asDouble());
      printf("SET: %20s -- %.3f\n","Q Value [MeV]",mdm->GetQValue());
    }
		else if(it.key().asString() == "excitationEnergy") {
      mdm->SetResidualEnergy(it->asDouble());
      printf("SET: %20s -- %.3f\n","Excitation Energy [MeV]",mdm->GetResidualEnergy());
		}
		else if(it.key().asString() == "excitationIsScattered") {
			exIsScattered = it->asBool();
			std::string sss = exIsScattered ? "YES":"NO";
			printf("SET: Excitation energy is of scattered particle? -- %s\n",sss.c_str());
		}
		// else if(it.key().asString() == "scatteredEnergy") { // ION ENERGY [MeV]
		// 	mdm->SetScatteredEnergy(it->asDouble());
		// 	printf("SET: %20s -- %.3f\n","Scattered Energy [MeV]",mdm->GetScatteredEnergy());
		// }
		else if(it.key().asString() == "scatteredAngles") { // ION ANGLES [deg]
			if(it->size() != 2) {
				std::cerr << "ERROR: Need two and only two scattering angles (low and high).\n";
				return 1;
			}
      for(unsigned int i = 0;i<it->size();i++) {
				scatteredAngles.push_back((*it)[i].asDouble());
      }
      printf("SET: %20s -- ","Scattered Angles [deg]");
      for(unsigned int i = 0;i<scatteredAngles.size();i++) {
				printf("%.3f ",scatteredAngles[i]);
      } printf("\n");
		}
		else if(it.key().asString() == "numberOfRays") { // Number of rays to send
			numberOfRays = it->asInt();
			printf("SET: %20s -- %li\n", "Number of Rays", numberOfRays);
		}
#if 0
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
#endif
		else if(it.key().asString() == "outputFile") {
			printf("Saving output to file: \"%s\"\n", it->asString().c_str());
			outputFileName = it->asString().c_str();
		}
		else if(it.key().asString() != "rayinFile") {
			fprintf(stderr, "WARNING: skipping unrecognized JSON field: %s\n",
							it.key().asString().c_str());
		}
	}
	// Set field values
	if(dipoleField.size() != 1) {
		std::cerr << "ERROR: Dipole field not set (need 1 and only 1 field)!\n";
		return 1;
	}
	if(multipoleField.empty()) {
		std::cout << "Using standard scaling for entrance multipole field...\n";
		mdm->SetMDMDipoleField(dipoleField.at(0));
	}
	else if(multipoleField.size() != dipoleField.size()) {
		std::cerr << "ERROR: Need the same number of multipole field values "
							<< "as dipole field values!\n";
		return 1;
	}
	else {
		std::cout << "Using manual value for entrance multipole field...\n";
		mdm->SetMDMDipoleMultipoleField(
			dipoleField.at(0), multipoleField.at(0));
	}
	
	// generate random rays
	auto do_ray_generation =
		[&]() {
			const double ThetaLow = scatteredAngles.at(0);
			const double ThetaHigh = scatteredAngles.at(1);
			const double deg = TMath::Pi()/180;
			const double amu = 931.49406;
			std::mt19937_64 rng;
			std::uniform_real_distribution<double> thetaDistribution(ThetaLow,ThetaHigh);
			std::uniform_real_distribution<double> phiDistribution(0,2*TMath::Pi());

			TFile fileout(outputFileName.c_str(),"recreate");
			TTree *tout = new TTree("mdm","MDM rays");
			double a1,x1,x2,x3,x4,b1,y1,y2,y3,y4,thcm,phcm,ekin,thlab,phlab,thx,thy;
			std::array<double*,17> addr =
				{&a1,&x1,&x2,&x3,&x4,&b1,&y1,&y2,&y3,&y4,
				 &thcm,&phcm,&ekin,&thlab,&phlab,&thx,&thy};
			std::array<std::string,17> nms =
				{"a1","x1","x2","x3","x4","b1","y1","y2","y3","y4",
				 "thcm","phcm","ekin","thlab","phlab","thx","thy"};
			for(int i=0; i< addr.size(); ++i){
				tout->Branch(nms.at(i).c_str(),addr.at(i),Form("%s/D",nms.at(i).c_str()));
			}
			int sent;tout->Branch("sent",&sent,"sent/I");
			int measured;tout->Branch("measured",&measured,"measured/I");
			TLorentzVector *lv3 = new TLorentzVector(0,0,0,0);
//			tout->Branch("lv3",&lv3);
			const double Q   = mdm->GetQValue();
			const double T1  = mdm->GetBeamEnergy();
			const double M1  = mdm->GetProjectileMass()*amu;
			const double P1  = sqrt(pow(M1+T1,2) - M1*M1);
			const double M2  = mdm->GetTargetMass()*amu;
			const double M3  = mdm->GetScatteredMass()*amu;
			const double Ex3 = exIsScattered ? mdm->GetResidualEnergy() : 0;
			const double Ex4 = exIsScattered ? 0 : mdm->GetResidualEnergy();
			const double M4 = M1 + M2 - (M3+Ex3) - Q - Ex4;

			TLorentzVector lv1(0,0,P1,M1+T1);
			TLorentzVector lv2(0,0,0,M2);
			TLorentzVector lv0 = lv1+lv2;
			const double S = lv0.Mag2();
			const double P3cm =
				sqrt((pow(S - pow(M3+Ex3,2) - pow(M4+Ex4,2), 2) -
							4*pow(M3+Ex3,2)*pow(M4+Ex4,2)) / (4*S));
			TVector3 bv = lv0.BoostVector();
			
			std::cout << "\nRunning Monte Carlo...\nNumber of rays left: ";
			std::flush(std::cout);
			while(numberOfRays-- > 0) {
				if((numberOfRays+1) % 100 == 0) {
					std::cout << numberOfRays+1 << "... ";
					std::flush(std::cout);
				}
				sent = measured = 0;
				a1=x1=x2=x3=x4=b1=y1=y2=y3=y4=thcm=phcm=ekin=thlab=phlab=thx=thy-1000;
				const double theta = thetaDistribution(rng)*deg;
				const double phi = phiDistribution(rng);
				thcm = theta/deg;
				phcm = phi/deg;
				lv3->SetPxPyPzE(P3cm*sin(theta)*cos(phi),
												P3cm*sin(theta)*sin(phi),
												P3cm*cos(theta),
												sqrt(P3cm*P3cm + pow(M3+Ex3,2))
					);
				lv3->Boost(bv);
				double thetaX = atan(lv3->Px()/lv3->Pz())/deg;
				double thetaY = atan(lv3->Py()/lv3->Pz())/deg;
				ekin = lv3->E()-lv3->M();
				thlab = lv3->Theta()/deg;
				phlab = lv3->Phi()/deg;
				thx = thetaX;
				thy = thetaY;
				if(fabs(thetaX) < 2 && fabs(thetaY) < 2){
					sent = 1;
					mdm->SetScatteredEnergy(lv3->E() - lv3->M());
					mdm->SetScatteredAngle(thetaX,thetaY);
					mdm->SendRay();
					mdm->GetOxfordWirePositions(a1,x1,x2,x3,x4,b1,y1,y2,y3,y4);
					a1 = (a1/1e3)/deg;
					b1 = (b1/1e3)/deg;
					if(fabs(x1) < 15) { measured = 1; }
				}
				tout->Fill();
			}
			tout->SetAlias("deg","TMath::Pi()/180");
			tout->SetAlias("z1","1.0*0");
			tout->SetAlias("z2","1.0*16.5");
			tout->SetAlias("z3","1.0*(16.5+17.7)");
			tout->SetAlias("z4","1.0*(16.5+17.7+15.9)");

			tout->Write();
			fileout.Close();
			return 0;
		};
		
	return do_ray_generation();
}

