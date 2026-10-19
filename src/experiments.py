"""Reproducible electrostatic experiments; every plotted number is exported."""
from pathlib import Path
from dataclasses import replace, asdict
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from electrostatics import (Device, FiniteVolume, EPS0, TERMINALS,
                            manufactured, vertical_exact)
from fem_reference import FiniteElement

plt.rcParams.update({"font.size":10,"axes.grid":True,"grid.alpha":0.22,
                     "figure.dpi":120,"savefig.dpi":170})


def write_csv(path, rows):
    rows=list(rows)
    with Path(path).open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def save(fig, folder, name):
    fig.savefig(folder/(name+".png"),bbox_inches="tight")
    plt.close(fig)


def basic(folder, h=0.5):
    m=FiniteVolume(Device(),h)
    vgs, vd_low, vd_high=0.5,0.05,0.7
    fig,axs=plt.subplots(2,1,figsize=(9,6),sharex=True,constrained_layout=True)
    potentials=[]
    for ax,vd in zip(axs,[vd_low,vd_high]):
        u=m.solve(vgs,vd); potentials.append(u)
        im=ax.pcolormesh(m.xe,m.ze,u,cmap="viridis",shading="flat",vmin=.4,vmax=1.6)
        ax.axhline(-5,color="white",ls="--",lw=1); ax.axhline(5,color="white",ls="--",lw=1)
        ax.set_ylabel("Vertical z (nm)")
        ax.set_title(f"2D depleted double-gate device: Vg={vgs:.2f} V, Vd={vd:.2f} V")
        fig.colorbar(im,ax=ax,label="Potential (V)")
    axs[-1].set_xlabel("Source-to-drain x (nm)")
    save(fig,folder,"potential_maps")
    np.savez_compressed(folder/"potential_maps.npz",x_edges_nm=m.xe,z_edges_nm=m.ze,
                        potential_low_V=potentials[0],potential_high_V=potentials[1])

    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    profiles=[]
    for ax,L in zip(axs,[30,70]):
        model=FiniteVolume(Device(length_nm=L),h)
        for vd in [vd_low,vd_high]:
            x,psi=model.centerline(model.solve(vgs,vd),model.voltages(vgs,vd))
            ax.plot(x,model.device.contact_V-psi,label=f"Vd={vd:.2f} V")
            b=model.barrier(vgs,vd)
            ax.plot(b["x_barrier_nm"],b["barrier_eV"],"o")
            profiles.extend(dict(length_nm=L,drain_V=vd,x_nm=xi,Ec_minus_Ec_source_eV=ei)
                            for xi,ei in zip(x,model.device.contact_V-psi))
        ax.set(xlabel="x (nm)",ylabel="Ec - Ec(source) (eV)",title=f"Centerline barrier, L={L} nm")
        ax.legend()
    save(fig,folder,"barrier_profiles"); write_csv(folder/"barrier_profiles.csv",profiles)

    c=m.capacitance(); np.savetxt(folder/"capacitance_F_per_m.csv",c,delimiter=",",
                                  header=",".join(TERMINALS),comments="")
    fig,axs=plt.subplots(1,3,figsize=(12,3.6),constrained_layout=True)
    im=axs[0].imshow(c*1e12,cmap="coolwarm")
    axs[0].set_xticks(range(4),["S","D","Gb","Gt"]); axs[0].set_yticks(range(4),["S","D","Gb","Gt"])
    axs[0].set_title("Maxwell capacitance (pF/m)")
    for j in range(4):
        for i in range(4): axs[0].text(i,j,f"{c[j,i]*1e12:.1f}",ha="center",va="center",fontsize=8)
    fig.colorbar(im,ax=axs[0],shrink=.75)
    for ax,field,title in [(axs[1],m.response[:,1].reshape(m.nz,m.nx),"Drain influence: d psi / d Vd"),
                           (axs[2],(m.response[:,2]+m.response[:,3]).reshape(m.nz,m.nx),"Tied-gate influence: d psi / d Vg")]:
        im=ax.pcolormesh(m.xe,m.ze,field,vmin=0,vmax=1,cmap="magma",shading="flat")
        ax.set(xlabel="x (nm)",ylabel="z (nm)",title=title)
        fig.colorbar(im,ax=ax,shrink=.75)
    save(fig,folder,"capacitance_influence")
    u=potentials[1]; volt=m.voltages(vgs,vd_high)
    charges=m.terminal_charges(u,volt)
    fixed=float(np.sum(m.rho*m.volume_nm2)*1e-18)
    return {"device":asdict(m.device),"h_nm":h,"barrier_0p5V_0p7V":m.barrier(vgs,vd_high),
            "dibl":m.dibl(),"terminal_order":list(TERMINALS),
            "terminal_charge_C_per_m":charges.tolist(),"fixed_body_charge_C_per_m":fixed,
            "relative_gauss_error":float(abs(charges.sum()+fixed)/max(np.max(abs(charges)),abs(fixed))),
            "relative_capacitance_asymmetry":float(np.max(abs(c-c.T))/np.max(abs(c))),
            "relative_capacitance_row_sum":float(np.max(abs(c.sum(1)))/np.max(abs(c))),
            "min_capacitance_eigenvalue_F_per_m":float(np.linalg.eigvalsh(c)[0]),
            "normalized_linear_residual":m.residual(u,volt)}


def geometry_sweeps(folder,h=0.5):
    lengths=np.arange(20,101,5,dtype=float)
    rows=[]
    for tsi in [6,10,16]:
        for L in lengths:
            d=Device(length_nm=L,silicon_nm=tsi)
            rows.append(dict(length_nm=L,silicon_nm=tsi,oxide_nm=d.oxide_nm,
                             **FiniteVolume(d,h).dibl()))
    write_csv(folder/"length_sweep.csv",rows)
    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for tsi in [6,10,16]:
        selected=[r for r in rows if r["silicon_nm"]==tsi]
        x=[r["length_nm"] for r in selected]
        axs[0].semilogy(x,[r["dibl_mV_V"] for r in selected],"o-",label=f"Si={tsi} nm")
        axs[1].plot(x,[r["threshold_low_V"] for r in selected],label=f"Si={tsi} nm, Vd=.05")
        axs[1].plot(x,[r["threshold_high_V"] for r in selected],ls="--",label=f"Si={tsi} nm, Vd=.70")
    axs[0].axhline(50,color="gray",ls="--",label="Design goal: 50 mV/V")
    axs[0].set(xlabel="Channel length (nm)",ylabel="Barrier-defined DIBL (mV/V)",title="Thinner silicon improves gate control")
    axs[1].set(xlabel="Channel length (nm)",ylabel="Barrier-defined Vth (V)",title="Criterion: 0.20 eV centerline barrier")
    for ax in axs: ax.legend(fontsize=8)
    save(fig,folder,"short_channel_scaling")

    ls=np.arange(20,71,5,dtype=float); ts=np.array([6,8,10,12,14,16],dtype=float)
    grid=np.empty((len(ts),len(ls))); design=[]
    for j,tsi in enumerate(ts):
        for i,L in enumerate(ls):
            model=FiniteVolume(Device(length_nm=L,silicon_nm=tsi),h)
            result=model.dibl(); grid[j,i]=result["dibl_mV_V"]
            design.append(dict(length_nm=L,silicon_nm=tsi,**result))
    write_csv(folder/"design_sweep.csv",design)
    fig,ax=plt.subplots(figsize=(8,4.5),constrained_layout=True)
    im=ax.pcolormesh(ls,ts,grid,cmap="viridis_r",shading="nearest")
    cs=ax.contour(ls,ts,grid,levels=[50],colors="white",linewidths=2)
    ax.clabel(cs,fmt="50 mV/V",inline=True)
    ax.set(xlabel="Channel length (nm)",ylabel="Silicon thickness (nm)",title="Electrostatic design space; tox=2 nm, kox=3.9")
    fig.colorbar(im,ax=ax,label="DIBL (mV/V)")
    save(fig,folder,"design_space")

    searches=[]
    for spacing in [h,h/2]:
        for tsi in [6,10,16]:
            def objective(L):
                return FiniteVolume(Device(length_nm=L,silicon_nm=tsi),spacing).dibl()["dibl_mV_V"]-50
            root=brentq(objective,20,100,xtol=0.01)
            searches.append(dict(silicon_nm=tsi,h_nm=spacing,minimum_length_nm=root,
                                 dibl_mV_V=objective(root)+50,goal_mV_V=50.0))
    write_csv(folder/"design_search.csv",searches)
    return searches


def dielectric_sweep(folder,h=0.5):
    cases=[(2,3.9),(4,7.8),(10,19.5)]
    rows=[]
    for tox,k in cases:
        for L in [20,30,40,50,70,100]:
            m=FiniteVolume(Device(length_nm=L,oxide_nm=tox,oxide_k=k),h)
            c=m.capacitance()
            # Charge on BOTH gates when BOTH are changed by the same amount.
            tied=float(c[2:4,2:4].sum())
            rows.append(dict(length_nm=L,oxide_nm=tox,oxide_k=k,EOT_nm=tox*3.9/k,
                             tied_gate_capacitance_F_per_m=tied,**m.dibl()))
    write_csv(folder/"dielectric_sweep.csv",rows)
    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for tox,k in cases:
        r=[v for v in rows if v["oxide_nm"]==tox]
        label=f"tox={tox} nm, k={k:g}"
        axs[0].semilogy([v["length_nm"] for v in r],[v["dibl_mV_V"] for v in r],"o-",label=label)
        axs[1].plot([v["length_nm"] for v in r],np.array([v["tied_gate_capacitance_F_per_m"] for v in r])*1e12,"o-",label=label)
    axs[0].set(xlabel="Channel length (nm)",ylabel="DIBL (mV/V)",title="Same EOT=2 nm; different lateral fringing")
    axs[1].set(xlabel="Channel length (nm)",ylabel="Tied-gate capacitance (pF/m)",title="Depleted geometric capacitance only")
    for ax in axs: ax.legend(fontsize=8)
    save(fig,folder,"same_eot_comparison")


def screening(folder,h=0.5):
    m=FiniteVolume(Device(length_nm=120),h)
    dpsi=m.response[:,1].reshape(m.nz,m.nx)
    x,value=m.centerline(dpsi,np.array([0,1,0,0]))
    distance=120-x
    fitmask=(distance>=15)&(distance<=45)
    slope,intercept=np.polyfit(distance[fitmask],np.log(value[fitmask]),1)
    fitted=-1/slope
    write_csv(folder/"screening.csv",(dict(distance_from_drain_nm=d,influence=v,
        exponential_fit=np.exp(intercept+slope*d)) for d,v in zip(distance,value)))
    exact=m.device.eigen_length_nm; approx=m.device.scale_length_nm
    result={"fit_interval_nm":[15,45],"fitted_length_nm":float(fitted),
            "layered_eigenmode_length_nm":exact,"parabolic_estimate_nm":approx,
            "relative_fit_error":float(abs(fitted/exact-1)),
            "relative_parabolic_error":float(approx/exact-1)}
    fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    visible=(distance>=1)&(distance<=65)
    axs[0].semilogy(distance[visible],value[visible],label="2D FVM drain response")
    axs[0].semilogy(distance[visible],np.exp(intercept+slope*distance[visible]),"--",label="Fit at 15-45 nm")
    axs[0].axvspan(15,45,color="gray",alpha=.12)
    axs[0].set(xlabel="Distance from drain (nm)",ylabel="d psi(centerline) / d Vd",title="Drain perturbation decays into the channel")
    axs[0].legend(fontsize=8)
    bars=axs[1].bar(["2D fitted","Layered mode","Parabolic"],[fitted,exact,approx],color=["#247ba0","#70c1b3","#f3b562"])
    for bar,v in zip(bars,[fitted,exact,approx]): axs[1].text(bar.get_x()+bar.get_width()/2,v+.08,f"{v:.4f}",ha="center")
    axs[1].set(ylabel="Scale length (nm)",title="An approximation has a measurable error",ylim=(0,7.5))
    save(fig,folder,"screening_length")
    return result


def verify(folder,h=0.5):
    mfg=[]; reference=[]
    for spacing in [2*h,h,h/2,h/4]:
        m=FiniteVolume(Device(length_nm=40),spacing)
        xx,zz=np.meshgrid(m.x,m.z); exact,rho=manufactured(m.device,xx,zz)
        actual=m.solve_custom(rho,lambda t,x,z: manufactured(m.device,x,z)[0])
        mfg.append(dict(h_nm=spacing,max_error_V=float(np.max(abs(actual-exact))),
                        rms_error_V=float(np.sqrt(np.mean((actual-exact)**2)))))
        m=FiniteVolume(Device(),spacing); fem=FiniteElement(Device(),spacing)
        u=m.solve(.5,.7); other=fem.solve(.5,.7)
        difference=u-fem.interpolate(other,m.x,m.z)
        xx,zz=np.meshgrid(m.x,m.z)
        interior=(xx>3)&(xx<27)
        reference.append(dict(h_nm=spacing,FVM_barrier_eV=m.barrier(.5,.7)["barrier_eV"],
                              FEM_barrier_eV=fem.barrier(other)["barrier_eV"],
                              potential_rms_difference_V=float(np.sqrt(np.mean(difference**2))),
                              full_domain_max_difference_V=float(np.max(abs(difference))),
                              contact_excluded_max_difference_V=float(np.max(abs(difference[interior])))))
    write_csv(folder/"manufactured_convergence.csv",mfg)
    write_csv(folder/"independent_mesh.csv",reference)
    m=FiniteVolume(Device(length_nm=300),h)
    u=m.solve(.7,.05)[:,m.nx//2]; analytical=vertical_exact(m.device,m.z,.7)
    write_csv(folder/"long_channel.csv",(dict(z_nm=z,numeric_V=a,analytical_V=b) for z,a,b in zip(m.z,u,analytical)))
    fig,axs=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    hs=np.array([r["h_nm"] for r in mfg]); err=np.array([r["max_error_V"] for r in mfg])
    axs[0,0].loglog(hs,err,"o-",label="Layered manufactured solution")
    axs[0,0].loglog(hs,err[-1]*(hs/hs[-1])**2,"--",label="Second-order guide")
    axs[0,0].set(xlabel="Mesh spacing (nm)",ylabel="Max error (V)",title="Known exact 2D solution; dielectric jump")
    axs[0,0].legend(fontsize=8)
    for key,label in [("FVM_barrier_eV","Finite volumes"),("FEM_barrier_eV","Triangular P1 FEM")]:
        axs[0,1].plot(hs,[r[key] for r in reference],"o-",label=label)
    axs[0,1].set(xlabel="Mesh spacing (nm)",ylabel="Barrier (eV)",title="Independent discretizations approach one limit")
    axs[0,1].legend(fontsize=8)
    for key,label in [("full_domain_max_difference_V","Full-domain maximum"),("contact_excluded_max_difference_V","Exclude x<3 or x>27 nm"),("potential_rms_difference_V","Full-domain RMS")]:
        axs[1,0].loglog(hs,[r[key] for r in reference],"o-",label=label)
    axs[1,0].set(xlabel="Mesh spacing (nm)",ylabel="FVM - FEM potential difference (V)",title="Contact endpoints reduce regularity")
    axs[1,0].legend(fontsize=7)
    for ax in [axs[0,0],axs[1,0]]:
        ax.set_xticks(hs,labels=[f"{v:g}" for v in hs])
        ax.tick_params(axis="x",which="minor",bottom=False,labelbottom=False)
    axs[1,1].plot(m.z,u,"o",ms=3,label="2D mid-channel, L=300 nm")
    axs[1,1].plot(m.z,analytical,label="1D layered analytical limit")
    axs[1,1].set(xlabel="Vertical z (nm)",ylabel="Potential (V)",title="Long-channel limit")
    axs[1,1].legend(fontsize=8)
    save(fig,folder,"verification")
    order=np.log(err[-2]/err[-1])/np.log(hs[-2]/hs[-1])
    return {"manufactured_finest_max_error_V":float(err[-1]),"manufactured_final_order":float(order),
            "finest_reference":reference[-1],"long_channel_max_error_V":float(np.max(abs(u-analytical)))}


def run(output, h=0.5):
    out=Path(output).resolve(); out.mkdir(parents=True,exist_ok=True)
    summary=basic(out,h)
    summary["design_search"]=geometry_sweeps(out,h)
    dielectric_sweep(out,h)
    summary["screening"]=screening(out,h)
    validation=verify(out,h)
    validation.update({k:summary[k] for k in ["relative_gauss_error","relative_capacitance_asymmetry",
                       "relative_capacitance_row_sum","normalized_linear_residual"]})
    summary["scope"]="Fixed-charge fully depleted electrostatics; barrier-defined threshold; no drain-current model or experimental data."
    for name,obj in [("summary.json",summary),("validation.json",validation)]:
        (out/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"DIBL_mV_V":summary["dibl"]["dibl_mV_V"],
                      "screening":summary["screening"],"validation":validation},indent=2))
    print(f"Saved eight figures, CSV/NPZ data and JSON to {out}")
    return summary,validation
