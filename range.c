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

void range()
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
	
	float Ei = 600; //enegia inizale del fascio in MeV
	float Ef = 0; //enegia finale del fascio in MeV
	
	//Inizializzazione variabili che regolano la struttura del ciclo
	int n = 100; //numero intervalli
	float step = (Ei-Ef)/n; //ampiezza di ogni inttervallo
	float range = 0; //valore iniziale del range che si accumula durante il ciclo per restituire il valore complessivo
	float drange; //spazio infinitesimo percorso dalla particella nelll'attraversare lo spessore di materiale
	float E_mean, gam, beta, s, eta, W_max, SP_mean, sigma, C;
	float E_fascio = Ei;

	//Struttura del ciclo
	for (int i=0; i<n ; i++){
	    Ei = E_fascio-i*step;
	    Ef = Ei-step;
	    E_mean = (Ei+Ef)/2;
		gam = 1+E_mean/E0;
		beta = TMath::Sqrt(1-TMath::Power((1/gam),2));
		s = E0e/E0;
		eta = beta*gam;
		W_max = (2*E0e*TMath::Power(eta,2))/(1+2*s*TMath::Sqrt(1+TMath::Power(eta,2))+TMath::Power(s,2));
		sigma = density_correction(eta, X0, sigma0, X1, C0, a, m);
		C = shell_correction(eta, I);
		SP_mean = K*ro*(Z/A)*TMath::Power((z/beta),2)*(TMath::Log((2*E0e*TMath::Power(eta,2)*W_max)/TMath::Power(I*1e-6,2))-2*TMath::Power(beta,2)-sigma-2*(C/Z));
		drange = (Ei-Ef)/SP_mean;
		range = range+drange;
		
	}

	//Stampa del risultato
	cout << "Il range di protoni a 600 MeV nel rame è pari a: " << range << "cm." << endl;
	
}