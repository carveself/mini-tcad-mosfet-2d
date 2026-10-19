"""Physical and numerical checks, independent of stored result files."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from electrostatics import Device,FiniteVolume,EPS0,manufactured,vertical_exact
from fem_reference import FiniteElement


class ElectrostaticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=FiniteVolume()

    def test_gauss_charge_balance(self):
        m=self.model; u=m.solve(.5,.7)
        charge=m.terminal_charges(u,m.voltages(.5,.7))
        body=np.sum(m.rho*m.volume_nm2)*1e-18
        self.assertLess(abs(charge.sum()+body)/np.max(abs(charge)),1e-10)

    def test_capacitance_reciprocity_and_passivity(self):
        c=self.model.capacitance(); scale=np.max(abs(c))
        np.testing.assert_allclose(c,c.T,rtol=0,atol=scale*1e-11)
        self.assertGreater(np.linalg.eigvalsh(c)[0],-scale*1e-10)
        self.assertTrue(np.all(np.diag(c)>0))
        self.assertTrue(np.all(c[~np.eye(4,dtype=bool)]<=0))

    def test_gauge_invariance(self):
        m=self.model
        before=m.solve(.5,.7,vs=0,bottom_vg=.4)
        after=m.solve(.8,1.0,vs=.3,bottom_vg=.7)
        np.testing.assert_allclose(after,before+.3,rtol=0,atol=1e-12)
        np.testing.assert_allclose(m.capacitance().sum(1),0,rtol=0,atol=1e-22)

    def test_charge_derivative_matches_capacitance(self):
        m=self.model; step=1e-5
        v=m.voltages(.5,.7); u=m.solve(.5,.7)
        q=m.terminal_charges(u,v)
        changed=m.solve(.5,.7+step); v[1]+=step
        derivative=(m.terminal_charges(changed,v)-q)/step
        np.testing.assert_allclose(derivative,m.capacitance()[:,1],rtol=1e-7,atol=1e-20)

    def test_top_bottom_symmetry(self):
        u=self.model.solve(.5,.7)
        np.testing.assert_allclose(u,u[::-1],rtol=0,atol=2e-13)

    def test_source_drain_symmetry_at_zero_drain(self):
        u=self.model.solve(.5,0)
        np.testing.assert_allclose(u,u[:,::-1],rtol=0,atol=2e-13)

    def test_gate_and_drain_lower_barrier(self):
        m=self.model
        b=m.barrier(.5,.05)["barrier_eV"]
        self.assertLess(m.barrier(.6,.05)["barrier_eV"],b)
        self.assertLess(m.barrier(.5,.7)["barrier_eV"],b)
        self.assertLess(m.barrier(.5,.7)["x_barrier_nm"],15)

    def test_longer_and_thinner_improve_dibl(self):
        baseline=self.model.dibl()["dibl_mV_V"]
        longer=FiniteVolume(Device(length_nm=50)).dibl()["dibl_mV_V"]
        thinner=FiniteVolume(Device(silicon_nm=6)).dibl()["dibl_mV_V"]
        self.assertGreater(baseline,0)
        self.assertLess(longer,baseline); self.assertLess(thinner,baseline)

    def test_independent_finite_element_barrier(self):
        m=self.model; fem=FiniteElement()
        value=fem.barrier(fem.solve(.5,.7))["barrier_eV"]
        self.assertLess(abs(m.barrier(.5,.7)["barrier_eV"]-value),0.001)

    def test_manufactured_dielectric_problem_second_order(self):
        errors=[]
        for h in [1,.5,.25]:
            m=FiniteVolume(Device(length_nm=40),h)
            xx,zz=np.meshgrid(m.x,m.z); expected,rho=manufactured(m.device,xx,zz)
            actual=m.solve_custom(rho,lambda t,x,z:manufactured(m.device,x,z)[0])
            errors.append(np.max(abs(expected-actual)))
        order=np.log(errors[-2]/errors[-1])/np.log(2)
        self.assertGreater(order,1.9); self.assertLess(order,2.1)

    def test_long_channel_analytical_profile(self):
        m=FiniteVolume(Device(length_nm=240))
        actual=m.solve(.7,.05)[:,m.nx//2]
        expected=vertical_exact(m.device,m.z,.7)
        self.assertLess(np.max(abs(actual-expected)),6e-6)

    def test_layered_eigenmode_matches_far_field_decay(self):
        m=FiniteVolume(Device(length_nm=120))
        _,v=m.centerline(m.response[:,1].reshape(m.nz,m.nx),np.array([0,1,0,0]))
        distance=120-np.r_[0,m.x,120]
        mask=(distance>15)&(distance<45)
        length=-1/np.polyfit(distance[mask],np.log(v[mask]),1)[0]
        self.assertLess(abs(length/m.device.eigen_length_nm-1),0.001)

    def test_no_charge_maximum_principle(self):
        m=FiniteVolume(Device(acceptors_cm3=0))
        u=m.solve(.5,.7,bottom_vg=.4)
        self.assertGreaterEqual(u.min(),.4-1e-12)
        self.assertLessEqual(u.max(),1.6+1e-12)

    def test_input_and_bracket_failures(self):
        for kwargs in [{"oxide_nm":0},{"oxide_k":-1},{"length_nm":np.nan},{"acceptors_cm3":-1}]:
            with self.assertRaises(ValueError): Device(**kwargs)
        with self.assertRaises(ValueError): FiniteVolume(h_nm=0)
        with self.assertRaises(ValueError): self.model.threshold(.7,bracket=(3,4))
        with self.assertRaises(ValueError): self.model.dibl(.7,.05)
        with self.assertRaises(ValueError): self.model.solve(np.nan,.05)


if __name__ == "__main__": unittest.main()
