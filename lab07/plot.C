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
  std::cout<<"Vectors are "<<N<<"-dimensional"<<std::endl;
  double valY;
  
  TMultiGraph* mg = new TMultiGraph();
  int n = 0;
  int nCorrect = 0;
  std::cout<<"Compare Eigenvalues"<<std::endl;
  std::cout<<" Numerical  | Analytical | difference"<<std::endl;
  std::cout<<"_____________________________________"<<std::endl;
  while(fValues>>eigenvalue){
    double En = omega *(n + 0.5);
    double error = TMath::Abs(En - eigenvalue);
    if(error<0.0001) nCorrect++;
    
    if(n<5){
      printf("%12.6lf|%12.6lf|%12.6lf\n",eigenvalue, En, error);
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
  }
  std::cout<<"Calculated "<<nCorrect<<" eigenvalues that were closer than 0.001"<<std::endl;

  TCanvas* c = new TCanvas("c", "Canva", 10, 10, w, h);
  c->DrawFrame(0,0, 400, 1);
  mg->Draw("l");
  mg->GetHistogram()->SetTitle("First five eigenvectors;x [a.u]; #Psi(x)");
  c->SaveAs("eigenvectors.png");
   
}
