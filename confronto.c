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


//Realizzazione grafico per confrontare l'andamento delle due espressione delle
//due formule di Bethe-Bloch in funzione dell'energia
void confronto(){
	
	//Inizializzazione dei parametri (nel caso di protone che impatta su un bersaglio di rame)
	
	float Z = 29; //numero atomico Cu
	float A = 64; //numero di massa Cu
	float z = 1; //carica particelle alfa in unità di carica elementare
	float ro = 8.9; //densità del Cu in g/cm^3
	float K = 0.1535; // costante in MeVcm^2/g
	float E0e = 0.511; //costante (massa a riposo e-) MeV/c^2
	float I = 322; //potenziale medio di ionizzazione del bersaglio in eV
	float E0 = 938.2; //energia a riposo protone in MeV
	float sigma0 = 0.08;
	float X0 = -0.0254;
	float X1 = 3.28;
	float C0 = -4.42;
	float a = 0.1434;
	float m = 2.90;
	
	//Inizializzazione variabili che regolano la struttura del ciclo
	float Emin = 10.0; //Energia minima del fascio in MeV
	float Emax = 1e5; //Energia massima del fascio in MeV
	int nstep = 1e3; //numero intervalli
	float dEstep = (Emax-Emin)/nstep;
	float gam, beta, s, eta, W_max, sigma, C;
	float E[nstep], SP[nstep], SP_corr[nstep];
	
	//Struttura del ciclo
	
	for (int i=0; i<nstep; i++){
		E[i] = Emin+i*dEstep;
		//Calcolo delle grandezze relativistiche
		gam = 1+E[i]/E0;
		beta = TMath::Sqrt(1-TMath::Power((1/gam),2));
		//Calcolo dell'energia massima trasferita in una singola collisione
		s = E0e/E0;
		eta = beta*gam;
		W_max = (2*E0e*TMath::Power(eta,2))/(1+2*s*TMath::Sqrt(1+TMath::Power(eta,2))+TMath::Power(s,2));
		//Calcolo dei termini correttivi
		sigma = density_correction(eta, X0, sigma0, X1, C0, a, m);
		C = shell_correction(eta, I);
		//Calcolo del valore di stopping power usando la formula senza e con  correzzioni 
		SP[i] = K*ro*(Z/A)*TMath::Power((z/beta),2)*(TMath::Log((2*E0e*TMath::Power(eta,2)*W_max)/TMath::Power(I*1e-6,2))-2*TMath::Power(beta,2));
		SP_corr[i] = K*ro*(Z/A)*TMath::Power((z/beta),2)*(TMath::Log((2*E0e*TMath::Power(eta,2)*W_max)/TMath::Power(I*1e-6,2))-2*TMath::Power(beta,2)-sigma-2*(C/Z));
		//cout << SP[i] << "	" << SP_corr[i] << endl;
		//cout << sigma << endl;
		
	}
	
	TCanvas *c1 = new TCanvas("c1","",800,600);
	c1->SetLogy();
	c1->SetLogx();
	//c1->SetGridx();
	//c1->SetGridy();
	
	//Creazione dei grafici
	
	TGraph *gr_SP =new TGraph(nstep,E,SP);
	gr_SP->SetTitle("senza correzioni");
	gr_SP->SetLineColor(kRed);
	gr_SP->SetMarkerColor(kRed);
	gr_SP->SetMarkerStyle(1);
	//gr_SP->Draw();
	
	TGraph *gr_SP_corr = new TGraph(nstep,E,SP_corr);
	gr_SP_corr->SetTitle("con correzioni");
	gr_SP_corr->SetLineColor(kBlue);
	gr_SP_corr->SetMarkerColor(kBlue);
	gr_SP_corr->SetMarkerStyle(1);
	//gr_SP_corr->Draw();
	
	//Creazione grafico su cui rappresentare le due curve
	TMultiGraph *mg = new TMultiGraph();

	//Opzioni grafiche
	mg->SetTitle("""; Energia (MeV); dE/dx (MeV/cm)"); 
	mg->GetXaxis()->CenterTitle(true);	
	mg->GetYaxis()->CenterTitle(true);
	mg->GetYaxis()->SetRangeUser(10,1e3);
	mg->GetXaxis()->SetLimits(Emin,Emax);
	
	//Grafici da rappresentare
	mg->Add(gr_SP,"cp");
	mg->Add(gr_SP_corr,"cp"); 
	mg->Draw("a");
	
	//Legenda
	c1->BuildLegend(0.9,0.75,0.48,0.9);
	
}