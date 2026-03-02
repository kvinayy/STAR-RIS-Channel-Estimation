# STAR-RIS Channel Estimation and Sum-Rate Maximization

This repository presents a study and implementation of **channel estimation and performance optimization for STAR-RIS assisted wireless communication systems**, targeting future **6G intelligent radio environments**.

The project combines **tensor-based channel estimation (PARAFAC)** with **optimization-driven beamforming and STAR-RIS configuration** to maximize system throughput.

---

## Motivation

Future 6G wireless systems require:

- High spectral efficiency
- Energy-efficient communication
- Reliable coverage in blocked environments

Traditional solutions rely on large antenna arrays, resulting in high hardware cost and power consumption.

**Reconfigurable Intelligent Surfaces (RIS)** enable programmable propagation environments by intelligently controlling signal reflection and transmission.

**STAR-RIS (Simultaneous Transmission and Reflection RIS)** further extends RIS capability by supporting users on both sides of the surface, enabling full-space wireless coverage.

---

## Project Overview

The implemented workflow consists of two major stages:

### Channel Estimation using PARAFAC

STAR-RIS phase switching during pilot transmission naturally forms a **3-way tensor**:

(BS antennas × Pilot symbols × RIS phase patterns)

Using **PARAFAC (CP tensor decomposition)**:

- Joint estimation of:
  - BS → STAR-RIS channel (**G**)
  - STAR-RIS → reflection users (**Hr**)
  - STAR-RIS → transmission users (**Ht**)
- Exploits sparsity and low-rank properties of mmWave channels
- Requires significantly fewer pilot signals
- Achieves NMSE performance close to oracle least-squares estimation

---

### Sum-Rate Maximization using PDD + BCD Optimization

Estimated channels are used to optimize system performance under STAR-RIS hardware constraints.

Optimization jointly determines:

- Base station beamforming matrix **W**
- STAR-RIS transmission coefficients
- STAR-RIS reflection coefficients
- Phase and amplitude coupling constraints

Techniques used:

- WMMSE Transformation
- Penalty Dual Decomposition (PDD)
- Block Coordinate Descent (BCD)

The algorithm converges to a **KKT optimal solution** for throughput maximization.

---

## System Model

- Base Station with **M antennas**
- STAR-RIS with **N programmable elements**
- Users divided into:
  - Transmission-side users
  - Reflection-side users
- Direct BS–User link assumed blocked
- mmWave channels modeled using **Saleh–Valenzuela sparse multipath model**

---

## Overall Workflow

Pilot Transmission
↓
Tensor Formation
↓
PARAFAC Channel Estimation
↓
Channel Reconstruction
↓
PDD + BCD Optimization
↓
Optimal Beamforming & STAR-RIS Configuration
↓
Maximum Sum-Rate

---

## Key Results

### PARAFAC Channel Estimation
- Accurate recovery of cascaded channels
- NMSE improves with increasing SNR
- Performance close to oracle LS baseline
- Reduced pilot overhead

### STAR-RIS Optimization
- Stable convergence of optimization algorithm
- Sum-rate ≈ **77–80 bit/s/Hz at SNR = 20 dB**
- Coupled transmission–reflection constraints satisfied

---

## Repository Structure

codes/ → Simulation & implementation codes
REPORT.pptx → Project presentation
IMPLEMENTATION-GUIDE.pdf → Implementation details
REFERENCES.md → Research papers used
README.md → Project documentation

---

## Technologies & Concepts

- Intelligent Reflecting Surfaces (IRS)
- STAR-RIS Systems
- mmWave Communication
- Tensor Decomposition (PARAFAC)
- Channel Estimation
- Beamforming Optimization
- PDD Optimization
- BCD Algorithms
- 6G Wireless Communication

---

## References

All referenced research papers are listed in [REFERENCES.md](REFERENCES.md)

---

## Authors

**Vinay Kusumanchi**  

**Rajdeep Alapati**
