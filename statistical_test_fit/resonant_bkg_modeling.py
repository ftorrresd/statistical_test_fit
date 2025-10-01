from ROOT import (  # type: ignore
    RooArgSet,  # type: ignore
    RooFit,  # type: ignore
    RooDataSet,  # type: ignore
    TFile,  # type: ignore
    RooWorkspace,  # type: ignore
)
from .fastplot import fastplot
import ws_helper 


def resonant_background_modeling_Higgs(input_file):
  ################################
  # Resonant Background Modeling #
  ################################

  w = RooWorkspace("resonant_background_ws")

  ### boson model

  w.factory("RooDoubleCB::resonant_background_model_boson("
            "boson_mass[75, 200], "
            "mean_boson[125, 75, 200], "
            "sigma_boson[2, 0.5, 4],"
            "alpha1_upsilon[3, 0, 10],"
            "n1_upsilon[0.5, 0.1, 50],"
            "alpha2_upsilon[3, 0, 10]," 
            "n2_upsilon[0.5, 0.1, 50]"
            ")")

  # upsilon model
  w.factory("RooBernstein::resonant_background_model_upsilon("
            "upsilon_mass["+str(config.configurations['dimuon_mass_range']['low'])+", "+str(config.configurations['dimuon_mass_range']['high'])+"],"
            "{"
            "1,"
            "p1[5, 0, 10]"
            "}"
            ")")


  # 2D model
  w.factory("PROD::resonant_background_model(resonant_background_model_boson,resonant_background_model_upsilon)")


  w.factory("weight[-100,100]")

  # load data
  f = TFile.Open(input_file)
  data = RooDataSet(
                      "resonant_background_data", 
                      "resonant_background_data", 
                      RooArgSet(w.var("boson_mass"), w.var("upsilon_mass"), w.var("weight")), 
                      RooFit.Import(f.Events),
                      RooFit.WeightVar(w.var("weight")),
                      )	
  getattr(w, 'import')(data)




  # fit to data      
  fit_result = w.pdf("resonant_background_model").fitTo(data, RooFit.Save())

  # plot data the and the pdf
  print "\n\n--> Saving plot " 
  nBins = 60
  w.var("boson_mass").SetTitle(r"m_{#mu#mu#gamma}")
  w.var("boson_mass").setUnit(r"GeV")
  w.var("upsilon_mass").SetTitle(r"m_{#mu#mu}")
  w.var("upsilon_mass").setUnit(r"GeV")
  # w.pdf("resonant_background_model_boson_gauss").SetTitle(r"Gaussian Component")
  # w.pdf("resonant_background_model_boson_cb").SetTitle(r"CB Component")
  # w.pdf("signal_model_upsilon").SetTitle(r"CB Component")

  # plot boson mass fit
  fastplot( w.pdf("resonant_background_model"), w.data("resonant_background_data"), w.var("boson_mass"), 
                    "outputs/resonant_background_fit/resonant_background_fit_boson" + inner_file_name + ".pdf",
                    # components=[
                    #             (w.pdf("resonant_background_model_boson_cb"), 10), 
                    #             (w.pdf("resonant_background_model_boson_gauss"), 10), 
                    #             ],
                    nbins=nBins,
                    legend=[0.6, 0.6, 0.9, 0.92], #legend=[0.2, 0.6, 0.5, 0.92],
                    )

  # plot upsilon mass fit
  fastplot( w.pdf("resonant_background_model"), w.data("resonant_background_data"), w.var("upsilon_mass"), 
                    "outputs/resonant_background_fit/resonant_background_fit_upsilon" + inner_file_name + ".pdf",
                    # components=[
                    #             (upsilon_1S, 10), 
                    #             (upsilon_2S, 10), 
                    #             (upsilon_3S, 10), 
                    #             (background, 20)
                    #             ],
                    nbins=nBins,
                    legend=[0.2, 0.7, 0.53, 0.9], #legend=[0.6, 0.2, 0.93, 0.4],
                    )

  # setting all var as constants
  w = ws_helper.set_constant(w)


  print("\n\n--> Fit parameters " )
  boson_parmeters = {}
  boson_parmeters["mean_bosonH"] = fit_result.floatParsFinal().find("mean_boson").getValV() #mean_boson.getValV()
  boson_parmeters["sigma_bosonH"] = fit_result.floatParsFinal().find("sigma_boson").getValV() #sigma_boson.getValV()
  boson_parmeters["alpha1_upsilonH"] = fit_result.floatParsFinal().find("alpha1_upsilon").getValV() #alpha1_upsilon.getValV()
  boson_parmeters["n1_upsilonH"] = fit_result.floatParsFinal().find("n1_upsilon").getValV() #n1_upsilon.getValV()
  boson_parmeters["alpha2_upsilonH"] = fit_result.floatParsFinal().find("alpha2_upsilon").getValV() #alpha2_upsilon.getValV()
  boson_parmeters["n2_upsilonH"] = fit_result.floatParsFinal().find("n2_upsilon").getValV() #n2_upsilon.getValV()



  print("data.sumEntries(): ",data.sumEntries())
