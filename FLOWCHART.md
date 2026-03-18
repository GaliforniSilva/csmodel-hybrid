# CS-Hybrid Model — Simulation Flowchart

This document provides a high-level scientific overview flowchart of the **CS-Hybrid 1D Coastal Profile Evolution Model** in both Mermaid and PlantUML formats, suitable for use as a figure in a scientific paper.

---

## A) Mermaid Diagram

```mermaid
flowchart TD
    START([Start]) --> A1["Model Initialisation\nCoastal profile geometry\nForcing time series · Initial sediment volumes"]

    A1 --> LOOP{{"Time step loop"}}

    LOOP --> B1["Hydrodynamic Forcing\nWave runup · Total Water Level"]
    B1   --> B2["Wave Erosion\nErosion and overwash fluxes"]
    B2   --> B3["Aeolian Transport\nBeach-to-dune sediment flux"]
    B3   --> B4["Subtidal Bar Evolution\nEquilibrium-based volume change"]
    B4   --> B5["Morphodynamic Profile Update\nDune geometry · Hard-structure interaction"]
    B5   --> B6["Vegetation Dynamics\nTwo-species seasonal growth model"]
    B6   --> B7["Sediment Budget Update\nBeach volume · Shoreline position"]

    B7 --> C1{"Simulation\ncomplete?"}
    C1 -- "No"  --> LOOP
    C1 -- "Yes" --> D1["Write Output\nTime series of all model variables"]

    D1 --> STOP([End])

    classDef init    fill:#d4ecd4,stroke:#2e7d32,color:#000
    classDef loop    fill:#fff9c4,stroke:#f9a825,color:#000
    classDef hydro   fill:#b3e0f7,stroke:#0277bd,color:#000
    classDef aeolian fill:#ffe0b2,stroke:#e65100,color:#000
    classDef bar     fill:#e8d5f5,stroke:#6a1b9a,color:#000
    classDef morpho  fill:#fce4ec,stroke:#880e4f,color:#000
    classDef veg     fill:#c8e6c9,stroke:#1b5e20,color:#000
    classDef vol     fill:#e0f2f1,stroke:#004d40,color:#000
    classDef out     fill:#f5f5f5,stroke:#424242,color:#000

    class START,A1 init
    class LOOP loop
    class B1,B2 hydro
    class B3 aeolian
    class B4 bar
    class B5 morpho
    class B6 veg
    class B7 vol
    class C1,D1,STOP out
```

---

## B) PlantUML Diagram

```plantuml
@startuml CS-Hybrid Simulation Flowchart

skinparam backgroundColor #FFFFFF
skinparam activity {
  BackgroundColor #F5F5F5
  BorderColor     #424242
  FontSize        12
}
skinparam ArrowColor #37474F
skinparam roundCorner 8

title CS-Hybrid 1D Coastal Profile Evolution Model

start

#d4ecd4:Model Initialisation\nCoastal profile geometry · Forcing time series · Initial sediment volumes;

while (Simulation running?) is (next step)

  #b3e0f7:Hydrodynamic Forcing\nWave runup · Total Water Level;

  #b3e0f7:Wave Erosion\nErosion and overwash fluxes;

  #ffe0b2:Aeolian Transport\nBeach-to-dune sediment flux;

  #e8d5f5:Subtidal Bar Evolution\nEquilibrium-based volume change;

  #fce4ec:Morphodynamic Profile Update\nDune geometry · Hard-structure interaction;

  #c8e6c9:Vegetation Dynamics\nTwo-species seasonal growth model;

  #e0f2f1:Sediment Budget Update\nBeach volume · Shoreline position;

endwhile (Simulation complete)

#d4ecd4:Write Output\nTime series of all model variables;

stop

@enduml
```

---

## Legend

### Module Colour Key

| Colour | Module | Description |
|--------|--------|-------------|
| 🟢 Green | **Initialisation / Output** | Profile setup, forcing loading, result writing |
| 🟡 Yellow | **Time loop** | Main integration loop |
| 🔵 Blue | **Hydrodynamics** | Wave runup, total water level, erosion and overwash |
| 🟠 Orange | **Aeolian transport** | Wind-driven beach-to-dune sediment flux |
| 🟣 Purple | **Subtidal bar** | Equilibrium bar volume and exchange with beach |
| 🔴 Red | **Morphodynamics** | Cross-shore profile update, hard-structure interaction |
| 🌿 Dark green | **Vegetation** | Two-species seasonal growth and decay |
| 🩵 Teal | **Sediment budget** | Beach volume balance and shoreline position |

---

### Key Variables

| Symbol | Meaning |
|--------|---------|
| **Hs** | Significant wave height [m] |
| **Tp** | Peak wave period [s] |
| **SWL** | Still water level (tide + surge) [m] |
| **TWL** | Total water level = wave runup + SWL [m] |
| **SLR** | Sea-level rise rate [m/yr] |
| **ρ₁, ρ₂** | Vegetation cover fraction, pioneer & mature species [-] |

