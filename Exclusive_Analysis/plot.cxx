#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include "iostream"
#include "TGraph.h"
#include "TLine.h"
#include <fstream>

#include <sstream>
#include <string>
#include <vector>

using namespace std;

// Function to split a string into substrings based on a delimiter
vector<string> splitString(const string &s, char delimiter) {
    vector<string> tokens;
    stringstream ss(s);
    string token;
    while (getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to read a CSV file and store its values in separate vectors
void readCSVFile(const string &filename, vector<double> &values1, vector<double> &values2, vector<double> &values3) {
    // Open the CSV file
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open the file " << filename << endl;
        return;
    }

    // Read and process each line of the file
    string line;
    while (getline(file, line)) {
        // Split the line into tokens
        vector<string> tokens = splitString(line, ',');

        // Ensure the line has exactly three values
        if (tokens.size() != 3) {
            cerr << "Error: Each line in file " << filename << " should have exactly three values." << endl;
            continue;
        }

        // Convert the tokens to doubles and push them into the respective vectors
        try {
            double val1 = stod(tokens[0]);
            double val2 = stod(tokens[1]);
            double val3 = stod(tokens[2]);
            values1.push_back(val1);
            values2.push_back(val2);
            values3.push_back(val3);
        } catch (...) {
            cerr << "Error: Failed to convert values to doubles in file " << filename << endl;
            continue;
        }
    }

    // Close the file
    file.close();
}

int main() {

    for (int n = 1; n <= 10; n++)
    {
            // Vectors to store values from each file
        vector<double> file1_values1, file1_values2, file1_values3;
        vector<double> file2_values1, file2_values2, file2_values3;
        vector<double> file3_values1, file3_values2, file3_values3;
        vector<double> file4_values1, file4_values2, file4_values3;
        vector<double> file5_values1, file5_values2, file5_values3;

        // Read data from each file
        readCSVFile("goodTof" +to_string(n)+".csv", file1_values1, file1_values2, file1_values3);
        readCSVFile("otherTpc" +to_string(n)+".csv", file2_values1, file2_values2, file2_values3);
        readCSVFile("true" +to_string(n)+".csv", file3_values1, file3_values2, file3_values3);
        readCSVFile("otherTpc" +to_string(n)+".csv", file2_values1, file2_values2, file2_values3);
        readCSVFile("pro" +to_string(n)+".csv", file4_values1, file4_values2, file4_values3);
        readCSVFile("tof" +to_string(n)+".csv", file5_values1, file5_values2, file5_values3);

     


        TCanvas *canvas = new TCanvas("canvas", "Graph with Markers", 900, 300);

        TGraph *graph = new TGraph();
        Int_t nPoints = file2_values1.size();
        for (Int_t i = 0; i < nPoints; ++i) {
        
            graph->SetPoint(i, file2_values1[i], file2_values2[i]);
            graph->SetMarkerStyle(24); // Set marker style (20 is a small circle)
            graph->SetMarkerColorAlpha(30,1);

            graph->SetMarkerSize(2.4 + 0.1 * file2_values3[i]); // Marker size increases with index
        }
            graph->GetXaxis()->SetTitle("#eta");
            graph->GetXaxis()->SetRangeUser(-9,9);
        graph->GetYaxis()->SetTitle("#phi [rad]");
        graph->GetXaxis()->SetTitleSize(0.1);
        graph->GetYaxis()->SetTitleSize(0.1);
        graph->GetXaxis()->SetLabelSize(0.1);
        graph->GetYaxis()->SetLabelSize(0.1);
    graph->GetXaxis()->SetTitleOffset(0.55);
    graph->GetYaxis()->SetTitleOffset(0.3);
        canvas->SetMargin(0.07, 0.045, 0.16, 0.01);


        TGraph *graph2 = new TGraph();
        Int_t nPoints2= file1_values2.size();
    
        for (Int_t i = 0; i < nPoints2; ++i) {
        
            graph2->SetPoint(i, file1_values1[i], file1_values2[i]);
            graph2->SetMarkerStyle(63); // Set marker style (20 is a small circle)
            
            graph2->SetMarkerColorAlpha(kRed,.7);
       
            graph2->SetMarkerSize(2.4 + 0.1 * file1_values3[i]); // Marker size increases with index


        }

        TGraph *graph3= new TGraph();
          
        
        Int_t nPoints3= file3_values1.size();
        for (Int_t i = 0; i < nPoints3; ++i) {
            graph3->SetPoint(i, file3_values1[i], file3_values2[i]);
          
            graph3->SetMarkerStyle(20); // Set marker style (20 is a small circle)
            graph3->SetMarkerColorAlpha(kBlue,0.5);
            graph3->SetMarkerSize(2.4 + 0.1 * file3_values3[i]); // Marker size increases with index




        }


        TGraph *graph4= new TGraph();
        Int_t nPoints4= file4_values1.size();
        for (Int_t i = 0; i < nPoints4; ++i) {
            graph4->SetPoint(i, file4_values1[i], file4_values2[i]);
           
            graph4->SetMarkerStyle(29); // Set marker style (20 is a small circle)
            graph4->SetMarkerColorAlpha(kOrange,1.0);
            graph4->SetMarkerSize(2.4 + 0.1 * file4_values3[i]); // Marker size increases with index
        }

    




        graph->GetXaxis()->SetLimits(-9,9); // Set range for x-axis from 0 to 11 for graph 1
        graph->GetYaxis()->SetRangeUser(-M_PI-0.1,M_PI+0.1); // Set range for y-axis from 0 to 20 for graph 2

        // Draw the first graph
        graph->Draw("AP"); // "A" to draw the axes, "P" to draw the markers
    // Draw the second graph on top of the first one
        graph3->Draw("P SAME"); // graph3"SAME" to draw on the same canvas
     
        // Set the range for the x-axis of the second graph
        graph2->GetXaxis()->SetLimits(-17,17); // Set range for x-axis from 4 to 14 for graph 2
        graph2->GetYaxis()->SetRangeUser(-M_PI-0.1,M_PI+0.1); // Set range for y-axis from 0 to 20 for graph 2

        // Draw the second graph on top of the first one
        graph2->Draw("P SAME"); // "SAME" to draw on the same canvas
        
        graph3->GetXaxis()->SetLimits(-17,17); // Set range for x-axis from 4 to 14 for graph 2
        graph3->GetYaxis()->SetRangeUser(-M_PI-0.1,M_PI+0.1); // Set range for y-axis from 0 to 20 for graph 2

    
        graph4->GetXaxis()->SetLimits(-17,17); // Set range for x-axis from 4 to 14 for graph 2
        graph4->GetYaxis()->SetRangeUser(-M_PI-0.1,M_PI+0.1); // Set range for y-axis from 0 to 20 for graph 2

        // Draw the second graph on top of the first one
        graph4->Draw("P SAME"); // "SAME" to draw on the same canvas
        



        TLatex* latex = new TLatex();
        latex->SetTextSize(0.1);
        latex->DrawLatexNDC(0.11, 0.2, "STAR pp #sqrt{s} = 510 GeV");


        TLegend *legend = new TLegend(0.65, 0.5, 1.04, 0.98);

        legend->SetMargin(0.1);
        legend->AddEntry(graph4, "intact protons", "P"); 
        legend->AddEntry(graph3, "true pions", "P"); 
        legend->AddEntry(graph2, "matched TPC tracks", "P"); 


            TGraph *graph5= new TGraph();
            Int_t nPoints5= file5_values1.size();
            for (Int_t i = 0; i < nPoints5; ++i)
            {
                graph5->SetPoint(i, file5_values1[i], file5_values2[i]);
   
                graph5->SetMarkerStyle(24); // Set marker style (20 is a small circle)
                graph5->SetMarkerColorAlpha(kBlack,1.0);
                graph5->SetMarkerSize(2.4 + 0.1 * file5_values3[i]); // Marker size increases with index
            }
            graph5->GetXaxis()->SetLimits(-17,17); // Set range for x-axis from 4 to 14 for graph 2
            graph5->GetYaxis()->SetRangeUser(-M_PI-0.1,M_PI+0.1); // Set range for y-axis from 0 to 20 for graph 2

            // Draw the second graph on top of the first one
            graph5->Draw("P SAME"); // "SAME" to draw on the same canvas
           // legend->AddEntry(graph5, "other TPC tracks ", "P"); 
  
        legend->AddEntry(graph5, "TPC tracks with TOF", "P"); 
        legend->AddEntry(graph, "other TPC tracks", "P"); 
        legend->SetBorderSize(0);
        legend->SetFillColorAlpha(0,0.0);
        legend->Draw();

        cout << n << endl;
   string str =  "topol_evt_" +to_string(n)+".pdf";

    const char* charPtr = str.c_str();
        canvas->SaveAs(charPtr);

    }
    return 0;
}
