//Funzione per calcolare la correzione di densità
float density_correction(float eta, float X0, float sigma0, float X1, float C0, float a, float m){
	//Calcolo della correzione di densità
	float X = TMath::Log10(eta);
	float sigma;
	if (X<=X0){
		sigma = sigma0*TMath::Power(10,2*(X-X0));
	}
		else if (X>X0 && X<X1){
			sigma = 4.6052*X+C0+a*TMath::Power(X1-X,m);
		}
		else{
			sigma = 4.6052*X+C0;
		}
		
	return sigma;
}

//Funzione per calcolare la correzione di shell
float shell_correction(float eta, float I){
	//Calcolo della correzione di shell
	float C;
	//if (eta>0.1){
	C = (0.422377*TMath::Power(eta,-2)+0.0304043*TMath::Power(eta,-4)-0.00038106*TMath::Power(eta,-6))*1e-6*TMath::Power(I,2) + 
		(3.850190*TMath::Power(eta,-2)-0.1667989*TMath::Power(eta,-4)+0.00157955*TMath::Power(eta,-6))*1e-9*TMath::Power(I,3);
	//}
	//else C = 0;
	return C;
}

//Realizzazione della curva di Bragg per un fascio 
//di protoni da 600 MeV che attraversa un bersaglio di rame

void curva_di_bragg()
{
	//Inizializzazione dei parametri
	float Z = 29; //numero atomico Cu
	float A = 64; //numero di massa Cu
	float z = 1; //carica particelle alpha in unità di carica elementare
	float ro = 8.9; //densità del Cu in g/cm^3
	float E0 = 938.2; //energia a riposo protone in MeV
	float K = 0.1535; // MeVcm^2/g
	float E0e = 0.511; //MeV/c^2
	float I = 322; //MeV
	float sigma0 = 0.08;
	float X0 = -0.0254;
	float X1 = 3.28;
	float C0 = -4.42;
	float a = 0.1434;
	float m = 2.90;
	
	float E_iniziale = 600; //enegia inizale del fascio in MeV
	float E_finale = 0; //enegia finale del fascio in MeV
	
	TCanvas *c1 = new TCanvas("c1","Curva di Bragg",800,600);
	c1->SetLogy();
	//c1->SetGridx();
	//c1->SetGridy();
	
	//Creazione del grafico
	TGraph *gr = new TGraph(); 
	gr->SetTitle("""; Spessore percorso dalla particella (cm); Stopping power (MeV/cm)");
	gr->GetXaxis()->CenterTitle(true);	
	gr->GetYaxis()->CenterTitle(true);
	
	//Inizializzazione variabili che regolano la struttura del ciclo 
	float E = 600;
	float range = 0;
	float gam = 1+E/E0;
	float beta, s, eta, W_max, SP, dE, sigma, C;
	float INC = 0.1; //incremento
	
	//Struttura del ciclo
	for (int i=0; E<E0 ; i++){
		gam = 1+E/E0;
		beta = TMath::Sqrt(1-TMath::Power((1/gam),2));
		s = E0e/E0;
		eta = beta*gam;
		W_max = (2*E0e*TMath::Power(eta,2))/(1+2*s*TMath::Sqrt(1+TMath::Power(eta,2))+TMath::Power(s,2));
		sigma = density_correction(eta, X0, sigma0, X1, C0, a, m);
		C = shell_correction(eta, I);
		SP =  K*ro*(Z/A)*TMath::Power((z/beta),2)*(TMath::Log((2*E0e*TMath::Power(eta,2)*W_max)/TMath::Power(I*1e-6,2))-2*TMath::Power(beta,2)-sigma-2*(C/Z));
	    dE = SP*INC;
		gr->SetPoint(i,range, SP);
		range = range+INC;
		E = E-dE;
	}

	//Opzioni grafiche
	gr->SetMarkerColor(kBlack);
	gr->SetMarkerStyle(1);
    gr->Draw("");
	
}