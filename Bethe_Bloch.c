//Funzione che restituisce il valore di dE/dx della formulla di Bethe Bloch una volta datole
//il valore del fattore di Lorentz e la massa a riposo della particella proiettile in MeV/c^2 
float mass_stopping_power_mean(float gam, float E0 ) //gam=fattore di Lorentz della particella ; E0=massa a riposo particella proiettile in MeV/c^2
{
	//Inizializzazione dei parametri (nel caso di protone che impatta su un bersaglio di rame)
	float Z = 29; //numero atomico Cu
	float A = 64; //numero di massa Cu
	float z = 1; //carica particelle alfa in unità di carica elementare
	float ro = 8.9; //densità del Cu in g/cm^3
	float K = 0.1535; // costante in MeVcm^2/g
	float E0e = 0.511; //costante (massa a riposo e-) MeV/c^2
	float I = 0.000322; //potenziale medio di eccitazione del bersaglio in MeV
	
	//Calcolo della grandezza relativistica beta 
	float beta = TMath::Sqrt(1-TMath::Power((1/gam),2));
   
	//Calcolo del valore medio di dE/dX dalla formula di Bethe-Bloch
	float s = E0e/E0;
	float eta = beta*gam;
	float W_max = (2*E0e*TMath::Power(eta,2))/(1+2*s*TMath::Sqrt(1+TMath::Power(eta,2))+TMath::Power(s,2));
	float SP_mean = K*ro*(Z/A)*TMath::Power((z/beta),2)*(TMath::Log((2*E0e*TMath::Power(eta,2)*W_max)/TMath::Power(I,2))-2*TMath::Power(beta,2));
	
	return SP_mean;
}

void Bethe_Bloch()
{	
	//----------------------Metodo n1
	
	cout << "Primo metodo: Valor medio." << endl;
	
	float E0 = 938.2; //energia a riposo protone in MeV
	float Ei = 600; //enegia inizale del fascio in MeV
	//float Ef = 500; //enegia finale del fascio in MeV
	float Ef = 0; //enegia finale del fascio in MeV (quando si calcola il range della particella)
	
	//Calcolo dell'energia media del fascio
	float E_mean = (Ei+Ef)/2;
	
	//Calcolo della grandezza relativistica gamma
	float gam = 1+E_mean/E0;
	
	//Calcolo del range medio
	float SP_mean = mass_stopping_power_mean(gam, E0);
	float x_mean = (Ei-Ef)/SP_mean;

	//Stampa del risultato
	cout << "Il valore è pari a circa " << x_mean << " cm." << endl;

	//----------------------Metodo n2
	
	cout << "Secondo metodo: Integrazione numerica." << endl;
	
	//Inizializzazione variabili che regolano la struttura del ciclo
	int n = 10; //numero intervalli
	float step = (Ei-Ef)/n; //ampiezza di ogni inttervallo
	float x = 0; //valore iniziale del range che si accumula durante il ciclo per restituire il valore complessivo
	float dx; //spazio infinitesimo percorso dalla particella nelll'attraversare lo spessore di materiale
	float E = Ei; //valore di enegia iniziale del fascio in MeV necessario per iniziare i cicli
	
	//Per stampare la procedura del ciclo
	cout << " Range" << "  " << "  SP" << "  " << " delta_x" << endl; //gli spazi servono per stampare meglio i risultati
	
	//Struttura del ciclo, rappresenta ciò che accade all'energia della particella quando attraversa lo spessore di materiale
	for (int i=0; i<n; i++){
	    Ei = E-i*step;
	    Ef = Ei-step;
	    E_mean = (Ei+Ef)/2;
	    gam = 1+E_mean/E0;
	    SP_mean = mass_stopping_power_mean(gam, E0);
		dx = (Ei-Ef)*(1/SP_mean);
	    x = x+dx;
		
		//Per stampare la procedura del ciclo
		cout << Ei << "-" <<  Ef << " " << SP_mean  << " " << dx << endl;
		
	}
	
	//Stampa del risultato
	//cout << "Lo spessore necessario risulta essere pari a " << x << " cm\n" << endl;
	cout << "Il range di protoni a 600 MeV nel rame è pari a: " << x << "cm." << endl; //(quando si calcola il range della particella)
		
	//----------------------Metodo n3
	
	TCanvas *c1 = new TCanvas("c1","Terzo metodo: estrapolazione grafica",800,600);
	//c1->SetLogy();
	c1->SetGridx();
	c1->SetGridy();
	
	//Creazione del grafico
	TGraph* gr = new TGraph(); 
	gr->SetTitle("""; Spessore percorso dalla particella (cm); Energia del fascio (MeV)"); 
	gr->GetXaxis()->CenterTitle(true);	
	gr->GetYaxis()->CenterTitle(true);
	
	//Inizializzazione variabili che regolano la struttura del ciclo 
	float range = 0;
	float SP, dE;
	float INC = 0.1; //incremento
	
	//cout << E << endl;
	//Struttura del ciclo
	gam = 1+E/E0; //valore iniziale di gamma
	for (int i=0; E>0; i++){
		gam = 1+E/E0;
		SP = mass_stopping_power_mean(gam, E0);
	    dE = SP*INC;
		gr->SetPoint(i,range, E);
		range = range+INC;
		E = E-dE;
		//cout <<  E << endl;
	}
	
	//Opzioni grafiche
	gr->SetMarkerColor(kBlack);
	gr->SetMarkerStyle(1);
    gr->Draw();

}