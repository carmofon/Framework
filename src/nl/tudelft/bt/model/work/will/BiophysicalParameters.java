/**
 * 
 */
package nl.tudelft.bt.model.work.will;

/**
 * @author wkc
 *
 */
public class BiophysicalParameters {

	public float maximumCellRadius = 10; // [um]
	
	public float minimumCellRadius = 1; // [um]
	
	public float oxygenDiffusivity = 6.4583e6f; // [um^2/h] (Y. Kim et al., 2011)
	
	public float EGFDiffusivity = 1.8648e5f; // [um^2/h] (Thorne, Hrabetova 2004)
	
	public float CSF1Diffusivity = 5.76e5f; // [um^2/h] (Wyckoff et al., 2004)
	
	public float oxygenTissueConcentration = 6.7e-6f; // [mol/cm^3] (Alexander, Mathematical medicine and biology, 2005)
	
	public float EGFTissueConcentration = 4e-7f; // [g/L] (Brzezinski, Lewinski, 1998)
	
	public float CSF1TissueConcentration = 0;
	
	public float EMTCoefficient = 0.0036f; // [1/h] (Y. Kim et al., 2011)
	
	public float biomassDensity = 3000; // [g/L]
	
	// relative maximum cell growth rate
	public float uMax = 0.0058f; // [1/h] (Y. Kim et al., 2011)
	
	public float Ko = 8.3e-3f; // [mol/m^3] 
	
	public float amoeboidSpeed = 24; // [um/h] (Goswami 2005) speed may be slower for mesenchymal and collective migration!
	
	// half-maximum response for gradient sensing accuracy
	public float EGFGradSensitivity = 0.02f; // [-] (Roussos 2011)
	
	public float CSF1GradSensitivity = 0.03f; // [-] (Tranquillo et al., 1988)
	
	
}
