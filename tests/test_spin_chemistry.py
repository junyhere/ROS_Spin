import unittest
import numpy as np
from spin_chemistry import EncounterParameters,P_DOUBLET,P_QUARTET,coherent_circuit_validation,downstream_ros,initial_density,propagate_encounter

class SpinChemistryTests(unittest.TestCase):
    def test_projectors(self):
        self.assertTrue(np.allclose(P_DOUBLET@P_DOUBLET,P_DOUBLET))
        self.assertTrue(np.allclose(P_QUARTET@P_QUARTET,P_QUARTET))
        self.assertTrue(np.allclose(P_DOUBLET@P_QUARTET,0))
        self.assertTrue(np.allclose(P_DOUBLET+P_QUARTET,np.eye(6)))
        self.assertAlmostEqual(float(np.trace(P_DOUBLET).real),2)
        self.assertAlmostEqual(float(np.trace(P_QUARTET).real),4)

    def test_manifold_initial_states(self):
        self.assertAlmostEqual(float(np.trace(P_DOUBLET@initial_density("doublet")).real),1)
        self.assertAlmostEqual(float(np.trace(P_QUARTET@initial_density("quartet")).real),1)
        self.assertAlmostEqual(float(np.trace(P_DOUBLET@initial_density("unpolarized")).real),1/3)

    def test_probability_accounting(self):
        out=propagate_encounter(EncounterParameters(k_doublet_s=2e6,k_quartet_s=2e5,k_escape_s=1e6),20e-6,samples=2001)
        total=out["superoxide_yield"]+out["escape_yield"]+out["unresolved_probability"]
        self.assertLess(abs(total-1),2e-4)

    def test_h2o2_is_downstream(self):
        out=downstream_ros(1e-6,1e-3,2e5,2e9,1e-6)
        self.assertEqual(out["hydrogen_peroxide_m"][0],0)
        self.assertGreater(out["hydrogen_peroxide_m"][-1],0)

    def test_circuit_embedding(self):
        p=EncounterParameters(field_t=(0,0,1e-4),exchange_rad_s=2e6,dipolar_rad_s=1e5,effective_hyperfine_rad_s=(3e5,0,0))
        self.assertLess(coherent_circuit_validation(p,1e-7),1e-12)

if __name__=="__main__": unittest.main()
