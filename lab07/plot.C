int w = 600;
int h = 600;

void plot(double omega = 0.001){
  // files storing eigenvalues and corresponding vectors
  ifstream fValues("eigenvalues.txt");
  ifstream fVectors("eigenvectors.txt");

  double eigenvalue;

  int N=-1; // dimension of the eigen-vectors
  fVectors>>N;
  if(N < 0){
    std::cout<<"Error while reading dimension of the vectors"<<std::endl;
    exit(1);
  }
  double valY;
  
  TMultiGraph* mg = new TMultiGraph();
  int n = 0;
  std::cout<<"Compare Eigenvalues"<<std::endl;
  std::cout<<" Numerical  | Analytical | difference"<<std::endl;
  std::cout<<"_____________________________________"<<std::endl;
  while(fValues>>eigenvalue){
    if(n<5){
      TGraph* gr = new TGraph();
      for(int i=0; i<N; i++){
	fVectors>>valY;
	gr->SetPoint(i, i, valY);
      }
      gr->SetLineColor(n+1);
      gr->SetTitle(Form("#lambda=%lf;x [a.u.];#Psi(x)",eigenvalue));
      mg->Add(gr);
    }
    n++;
    double En = omega *(n + 0.5);
    printf("%lf5.6|%lf5.6|%lf5.6",eigenvalue, En, TMath::Abs(En - eigenvalue));
  }

  TCanvas* c = new TCanvas("c", "Canva", 10, 10, w, h);
  c->DrawFrame(0,0, 400, 10);
  mg->Draw("l");
  mg->GetHistogram()->SetTitle("First five eigenvectors;x [a.u]; #Psi(x)");
  c->SaveAs("eigenvectors.png");
   
}
