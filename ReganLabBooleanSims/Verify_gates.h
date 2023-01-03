//
//  Verify_gates.h
//  Boolean_code_for_IS
//
//  Created by Erzsó Ravasz Regan on 7/18/18.
//  Copyright © 2018 Regan_Group. All rights reserved.
//

#ifndef Verify_gates_h
#define Verify_gates_h

// these were generated from the list of nodes using the program Patterns.
    // The reg exp is:    (.*)
// The replacement is:    #define $1 A->s[A->BRN->MDAT->get_node_ID("$1")]
//                        #define $1_n A->BRN->MDAT->get_node_ID("$1")


#define GF A->s[A->BRN->MDAT->get_node_ID("GF")]
#define GF_n A->BRN->MDAT->get_node_ID("GF")
#define RTK A->s[A->BRN->MDAT->get_node_ID("RTK")]
#define RTK_n A->BRN->MDAT->get_node_ID("RTK")
#define Grb2 A->s[A->BRN->MDAT->get_node_ID("Grb2")]
#define Grb2_n A->BRN->MDAT->get_node_ID("Grb2")
#define Ras A->s[A->BRN->MDAT->get_node_ID("Ras")]
#define Ras_n A->BRN->MDAT->get_node_ID("Ras")
#define RAF A->s[A->BRN->MDAT->get_node_ID("RAF")]
#define RAF_n A->BRN->MDAT->get_node_ID("RAF")
#define mTORC2 A->s[A->BRN->MDAT->get_node_ID("mTORC2")]
#define mTORC2_n A->BRN->MDAT->get_node_ID("mTORC2")
#define PI3K A->s[A->BRN->MDAT->get_node_ID("PI3K")]
#define PI3K_n A->BRN->MDAT->get_node_ID("PI3K")
#define PIP3 A->s[A->BRN->MDAT->get_node_ID("PIP3")]
#define PIP3_n A->BRN->MDAT->get_node_ID("PIP3")
#define PDK1 A->s[A->BRN->MDAT->get_node_ID("PDK1")]
#define PDK1_n A->BRN->MDAT->get_node_ID("PDK1")
#define AKT_B A->s[A->BRN->MDAT->get_node_ID("AKT_B")]
#define AKT_B_n A->BRN->MDAT->get_node_ID("AKT_B")
#define p110_H A->s[A->BRN->MDAT->get_node_ID("p110_H")]
#define p110_H_n A->BRN->MDAT->get_node_ID("p110_H")
#define PI3K_H A->s[A->BRN->MDAT->get_node_ID("PI3K_H")]
#define PI3K_H_n A->BRN->MDAT->get_node_ID("PI3K_H")
#define AKT_H A->s[A->BRN->MDAT->get_node_ID("AKT_H")]
#define AKT_H_n A->BRN->MDAT->get_node_ID("AKT_H")
#define FoxO3 A->s[A->BRN->MDAT->get_node_ID("FoxO3")]
#define FoxO3_n A->BRN->MDAT->get_node_ID("FoxO3")
#define PLCgamma A->s[A->BRN->MDAT->get_node_ID("PLCgamma")]
#define PLCgamma_n A->BRN->MDAT->get_node_ID("PLCgamma")
#define NeddL4 A->s[A->BRN->MDAT->get_node_ID("NeddL4")]
#define NeddL4_n A->BRN->MDAT->get_node_ID("NeddL4")
#define FoxO1 A->s[A->BRN->MDAT->get_node_ID("FoxO1")]
#define FoxO1_n A->BRN->MDAT->get_node_ID("FoxO1")
#define p21_mRNA A->s[A->BRN->MDAT->get_node_ID("p21_mRNA")]
#define p21_mRNA_n A->BRN->MDAT->get_node_ID("p21_mRNA")
#define TSC2 A->s[A->BRN->MDAT->get_node_ID("TSC2")]
#define TSC2_n A->BRN->MDAT->get_node_ID("TSC2")
#define PRAS40 A->s[A->BRN->MDAT->get_node_ID("PRAS40")]
#define PRAS40_n A->BRN->MDAT->get_node_ID("PRAS40")
#define Rheb A->s[A->BRN->MDAT->get_node_ID("Rheb")]
#define Rheb_n A->BRN->MDAT->get_node_ID("Rheb")
#define mTORC1 A->s[A->BRN->MDAT->get_node_ID("mTORC1")]
#define mTORC1_n A->BRN->MDAT->get_node_ID("mTORC1")
#define S6K A->s[A->BRN->MDAT->get_node_ID("S6K")]
#define S6K_n A->BRN->MDAT->get_node_ID("S6K")
#define eIF4E A->s[A->BRN->MDAT->get_node_ID("eIF4E")]
#define eIF4E_n A->BRN->MDAT->get_node_ID("eIF4E")
#define GSK3 A->s[A->BRN->MDAT->get_node_ID("GSK3")]
#define GSK3_n A->BRN->MDAT->get_node_ID("GSK3")
#define p21 A->s[A->BRN->MDAT->get_node_ID("p21")]
#define p21_n A->BRN->MDAT->get_node_ID("p21")
#define pRB A->s[A->BRN->MDAT->get_node_ID("pRB")]
#define pRB_n A->BRN->MDAT->get_node_ID("pRB")
#define p27Kip1 A->s[A->BRN->MDAT->get_node_ID("p27Kip1")]
#define p27Kip1_n A->BRN->MDAT->get_node_ID("p27Kip1")
#define Myc A->s[A->BRN->MDAT->get_node_ID("Myc")]
#define Myc_n A->BRN->MDAT->get_node_ID("Myc")
#define CyclinD1 A->s[A->BRN->MDAT->get_node_ID("CyclinD1")]
#define CyclinD1_n A->BRN->MDAT->get_node_ID("CyclinD1")
#define E2F1 A->s[A->BRN->MDAT->get_node_ID("E2F1")]
#define E2F1_n A->BRN->MDAT->get_node_ID("E2F1")
#define CyclinE A->s[A->BRN->MDAT->get_node_ID("CyclinE")]
#define CyclinE_n A->BRN->MDAT->get_node_ID("CyclinE")
#define ORC A->s[A->BRN->MDAT->get_node_ID("ORC")]
#define ORC_n A->BRN->MDAT->get_node_ID("ORC")
#define Cdc6 A->s[A->BRN->MDAT->get_node_ID("Cdc6")]
#define Cdc6_n A->BRN->MDAT->get_node_ID("Cdc6")
#define Cdt1 A->s[A->BRN->MDAT->get_node_ID("Cdt1")]
#define Cdt1_n A->BRN->MDAT->get_node_ID("Cdt1")
#define Pre_RC A->s[A->BRN->MDAT->get_node_ID("Pre-RC")]
#define Pre_RC_n A->BRN->MDAT->get_node_ID("Pre-RC")
#define geminin A->s[A->BRN->MDAT->get_node_ID("geminin")]
#define geminin_n A->BRN->MDAT->get_node_ID("geminin")
#define CyclinA_mRNA A->s[A->BRN->MDAT->get_node_ID("CyclinA_mRNA")]
#define CyclinA_mRNA_n A->BRN->MDAT->get_node_ID("CyclinA_mRNA")
#define Emi1 A->s[A->BRN->MDAT->get_node_ID("Emi1")]
#define Emi1_n A->BRN->MDAT->get_node_ID("Emi1")
#define FoxM1 A->s[A->BRN->MDAT->get_node_ID("FoxM1")]
#define FoxM1_n A->BRN->MDAT->get_node_ID("FoxM1")
#define Cdc25A A->s[A->BRN->MDAT->get_node_ID("Cdc25A")]
#define Cdc25A_n A->BRN->MDAT->get_node_ID("Cdc25A")
#define CyclinA A->s[A->BRN->MDAT->get_node_ID("CyclinA")]
#define CyclinA_n A->BRN->MDAT->get_node_ID("CyclinA")
#define Wee1 A->s[A->BRN->MDAT->get_node_ID("Wee1")]
#define Wee1_n A->BRN->MDAT->get_node_ID("Wee1")
#define UbcH10 A->s[A->BRN->MDAT->get_node_ID("UbcH10")]
#define UbcH10_n A->BRN->MDAT->get_node_ID("UbcH10")
#define CyclinB A->s[A->BRN->MDAT->get_node_ID("CyclinB")]
#define CyclinB_n A->BRN->MDAT->get_node_ID("CyclinB")
#define Cdc25B A->s[A->BRN->MDAT->get_node_ID("Cdc25B")]
#define Cdc25B_n A->BRN->MDAT->get_node_ID("Cdc25B")
#define Plk1 A->s[A->BRN->MDAT->get_node_ID("Plk1")]
#define Plk1_n A->BRN->MDAT->get_node_ID("Plk1")
#define Cdc25C A->s[A->BRN->MDAT->get_node_ID("Cdc25C")]
#define Cdc25C_n A->BRN->MDAT->get_node_ID("Cdc25C")
#define Cdk1 A->s[A->BRN->MDAT->get_node_ID("Cdk1")]
#define Cdk1_n A->BRN->MDAT->get_node_ID("Cdk1")
#define pAPC A->s[A->BRN->MDAT->get_node_ID("pAPC")]
#define pAPC_n A->BRN->MDAT->get_node_ID("pAPC")
#define Cdc20 A->s[A->BRN->MDAT->get_node_ID("Cdc20")]
#define Cdc20_n A->BRN->MDAT->get_node_ID("Cdc20")
#define Cdh1 A->s[A->BRN->MDAT->get_node_ID("Cdh1")]
#define Cdh1_n A->BRN->MDAT->get_node_ID("Cdh1")
#define Replication A->s[A->BRN->MDAT->get_node_ID("Replication")]
#define Replication_n A->BRN->MDAT->get_node_ID("Replication")
#define f4N_DNA A->s[A->BRN->MDAT->get_node_ID("4N_DNA")]
#define f4N_DNA_n A->BRN->MDAT->get_node_ID("4N_DNA")
#define U_Kinetochores A->s[A->BRN->MDAT->get_node_ID("U_Kinetochores")]
#define U_Kinetochores_n A->BRN->MDAT->get_node_ID("U_Kinetochores")
#define Mad2 A->s[A->BRN->MDAT->get_node_ID("Mad2")]
#define Mad2_n A->BRN->MDAT->get_node_ID("Mad2")
#define A_Kinetochores A->s[A->BRN->MDAT->get_node_ID("A_Kinetochores")]
#define A_Kinetochores_n A->BRN->MDAT->get_node_ID("A_Kinetochores")
#define Plk1_H A->s[A->BRN->MDAT->get_node_ID("Plk1_H")]
#define Plk1_H_n A->BRN->MDAT->get_node_ID("Plk1_H")
#define Ect2 A->s[A->BRN->MDAT->get_node_ID("Ect2")]
#define Ect2_n A->BRN->MDAT->get_node_ID("Ect2")
#define Casp8 A->s[A->BRN->MDAT->get_node_ID("Casp8")]
#define Casp8_n A->BRN->MDAT->get_node_ID("Casp8")
#define Casp2 A->s[A->BRN->MDAT->get_node_ID("Casp2")]
#define Casp2_n A->BRN->MDAT->get_node_ID("Casp2")
#define MCL_1 A->s[A->BRN->MDAT->get_node_ID("MCL-1")]
#define MCL_1_n A->BRN->MDAT->get_node_ID("MCL-1")
#define BCLXL A->s[A->BRN->MDAT->get_node_ID("BCLXL")]
#define BCLXL_n A->BRN->MDAT->get_node_ID("BCLXL")
#define BCL2 A->s[A->BRN->MDAT->get_node_ID("BCL2")]
#define BCL2_n A->BRN->MDAT->get_node_ID("BCL2")
#define BAD A->s[A->BRN->MDAT->get_node_ID("BAD")]
#define BAD_n A->BRN->MDAT->get_node_ID("BAD")
#define BIK A->s[A->BRN->MDAT->get_node_ID("BIK")]
#define BIK_n A->BRN->MDAT->get_node_ID("BIK")
#define BIM A->s[A->BRN->MDAT->get_node_ID("BIM")]
#define BIM_n A->BRN->MDAT->get_node_ID("BIM")
#define BID A->s[A->BRN->MDAT->get_node_ID("BID")]
#define BID_n A->BRN->MDAT->get_node_ID("BID")
#define BAK A->s[A->BRN->MDAT->get_node_ID("BAK")]
#define BAK_n A->BRN->MDAT->get_node_ID("BAK")
#define BAX A->s[A->BRN->MDAT->get_node_ID("BAX")]
#define BAX_n A->BRN->MDAT->get_node_ID("BAX")
#define Cyto_C A->s[A->BRN->MDAT->get_node_ID("Cyto_C")]
#define Cyto_C_n A->BRN->MDAT->get_node_ID("Cyto_C")
#define SMAC A->s[A->BRN->MDAT->get_node_ID("SMAC")]
#define SMAC_n A->BRN->MDAT->get_node_ID("SMAC")
#define IAPs A->s[A->BRN->MDAT->get_node_ID("IAPs")]
#define IAPs_n A->BRN->MDAT->get_node_ID("IAPs")
#define Casp9 A->s[A->BRN->MDAT->get_node_ID("Casp9")]
#define Casp9_n A->BRN->MDAT->get_node_ID("Casp9")
#define Casp3 A->s[A->BRN->MDAT->get_node_ID("Casp3")]
#define Casp3_n A->BRN->MDAT->get_node_ID("Casp3")
#define CAD A->s[A->BRN->MDAT->get_node_ID("CAD")]
#define CAD_n A->BRN->MDAT->get_node_ID("CAD")


#define GF_High A->s[A->BRN->MDAT->get_node_ID("GF_High")]
#define GF_High_n A->BRN->MDAT->get_node_ID("GF_High")
#define SOS A->s[A->BRN->MDAT->get_node_ID("SOS")]
#define SOS_n A->BRN->MDAT->get_node_ID("SOS")

#define ERK A->s[A->BRN->MDAT->get_node_ID("ERK")]
#define ERK_n A->BRN->MDAT->get_node_ID("ERK")
#define Ca2 A->s[A->BRN->MDAT->get_node_ID("Ca2+")]
#define Ca2_n A->BRN->MDAT->get_node_ID("Ca2+")
#define IP3 A->s[A->BRN->MDAT->get_node_ID("IP3")]
#define IP3_n A->BRN->MDAT->get_node_ID("IP3")
#define DAG A->s[A->BRN->MDAT->get_node_ID("DAG")]
#define DAG_n A->BRN->MDAT->get_node_ID("DAG")
#define CHK1 A->s[A->BRN->MDAT->get_node_ID("CHK1")]
#define CHK1_n A->BRN->MDAT->get_node_ID("CHK1")
#define DR4_5 A->s[A->BRN->MDAT->get_node_ID("DR4_5")]
#define DR4_5_n A->BRN->MDAT->get_node_ID("DR4_5")


//unsigned short int GF_function(Boolean_Dynamics *A){
//    return(a);
//}

unsigned int *ge=NULL;

void evaluate_gate_expressions(Boolean_Dynamics_TRACKER *A){

    ge[GF_n] = GF || GF_High ;
    ge[RTK_n] = !CAD && (GF_High || GF) ;
    ge[Grb2_n] = RTK && GF_High ;
    ge[Ras_n] = Grb2 && SOS ;
    ge[RAF_n] = !Casp3 && Ras ;
    ge[mTORC2_n] = PIP3 || !S6K ;
    ge[PI3K_n] = Ras || RTK ;
    ge[PIP3_n] = PI3K_H || PI3K ;
    ge[PDK1_n] = PI3K && PIP3 ;
    ge[AKT_B_n] = !Casp3 && PIP3 && (PDK1 || mTORC2) ;
    ge[p110_H_n] = (FoxO3 && !NeddL4) || (p110_H && (FoxO3 || !NeddL4)) ;
    ge[PI3K_H_n] = p110_H && RTK && PI3K && Ras ;
    ge[AKT_H_n] = AKT_B && p110_H && PI3K_H && Ras && PIP3 && PDK1 && mTORC2 ;
    ge[FoxO3_n] = !(AKT_B || AKT_H || ERK) || (  !(AKT_H && (Plk1 || Plk1_H || AKT_B || ERK)) && !(Plk1 && Plk1_H && ERK)) ;
    ge[PLCgamma_n] = RTK && Grb2 && p110_H && PI3K_H && PIP3 ;
    ge[NeddL4_n] = Ca2 && IP3 ;
    ge[FoxO1_n] = !Plk1 && !AKT_H ;
    ge[p21_mRNA_n] = (FoxO1 && FoxO3) || (!Myc && (FoxO1 || FoxO3)) ;
    ge[TSC2_n] = !AKT_H || !(AKT_B || ERK) ;
    ge[PRAS40_n] = !AKT_H && (!mTORC1 || !AKT_B) ;
    ge[Rheb_n] = !TSC2 && DAG ;
    ge[mTORC1_n] = !Casp3 && ( (Rheb && !PRAS40) || E2F1 || (CyclinB && Cdk1 && GSK3) ) ;
    ge[S6K_n] = !Casp3 && mTORC1 ;
    ge[eIF4E_n] = !Casp3 && mTORC1 ;
    ge[GSK3_n] =  !AKT_H && !(S6K && ERK) ;
    ge[p21_n] = p21_mRNA && !Casp3 && !CyclinE ;
    ge[pRB_n] = !Casp3 && !CyclinD1 && !CyclinA && (p27Kip1 || !CyclinE) ;
    ge[p27Kip1_n] = !Casp3 && !CyclinD1 && !(Cdk1 && CyclinB) && ( (!CyclinE && (FoxO3 && FoxO1)) || (!CyclinA && (FoxO3 || FoxO1) ) || (!CyclinE && !CyclinA)) ;
    ge[Myc_n] = (ERK && (eIF4E || !GSK3)) || ( E2F1 && !pRB && (eIF4E || ERK || !GSK3) ) ;
    ge[CyclinD1_n] = !CHK1 && ( ( !p21 && ( (!GSK3 && (Myc || E2F1)) || (Myc && CyclinD1) || (Myc && E2F1) || (E2F1 && CyclinD1) ) ) || ( !pRB && E2F1 && ( (Myc && CyclinD1) || (Myc && !GSK3) || (CyclinD1 && !GSK3) ) ) ) ;
    ge[E2F1_n] = !(CAD || CyclinA || pRB) && (E2F1 || Myc) ;
    ge[CyclinE_n] = E2F1 && Cdc6 && Pre_RC && !(pRB || p27Kip1 || CHK1 || Casp3) ;
    ge[ORC_n] = E2F1 || (Pre_RC && Cdt1 && Cdc6) ;
    ge[Cdc6_n] = !Casp3 && !(f4N_DNA && CyclinA) && ( (E2F1 && ORC) || (Pre_RC && ORC && Cdc6 && Cdt1) ) ;
    ge[Cdt1_n] = !geminin && ORC && Cdc6 && !(CyclinE && CyclinA && Cdc25A) && ( (Pre_RC && (E2F1 || Myc)) || (E2F1 && (Myc || !pRB) && (FoxO3 || !f4N_DNA))) ;
    ge[Pre_RC_n] = ORC && Cdc6 && Cdt1 && !(Replication && f4N_DNA) ;
    ge[geminin_n] = E2F1 && !Cdh1 && !(pAPC && Cdc20) ;
    ge[CyclinA_mRNA_n] = !CAD && ( (E2F1 && !pRB) || FoxM1 ) ;
    ge[Emi1_n] = (E2F1 || !pRB || !p21) && !(Plk1 && CyclinB && Cdk1 && (U_Kinetochores || A_Kinetochores)) ;
    ge[FoxM1_n] = (Myc && CyclinE) || (CyclinA && Cdc25A && Cdc25B) || (Plk1 && CyclinB && Cdk1) ;
    ge[Cdc25A_n] = ( (FoxM1 && E2F1 && !pRB) || (!Cdh1 && (FoxM1 || (E2F1 && !pRB) ) ) ) && (!(GSK3 || CHK1) || CyclinE || CyclinA || (CyclinB && Cdk1) )  ;
                    
    ge[CyclinA_n] = CyclinA_mRNA && !pAPC && ( (Cdc25A && (!Cdh1 || Emi1) ) || (CyclinA && ( (!Cdh1 && (Emi1 || !UbcH10) ) || (Emi1 && !UbcH10) ) ) ) ;
    ge[Wee1_n] = !Casp3 && !(Cdk1 && CyclinB) && (Replication || CHK1) && (CHK1 || !(Cdk1 && CyclinA && Plk1)) ;
    ge[UbcH10_n] = !Cdh1 || (UbcH10 && (Cdc20 || CyclinA || CyclinB) ) ;
    ge[CyclinB_n] = (FoxM1 || (FoxO3 && CyclinB) ) && !(Cdh1 || (pAPC && Cdc20)) ;
    ge[Cdc25B_n] = FoxM1 && f4N_DNA ;
    ge[Plk1_n] = !Cdh1 && (FoxM1 || Plk1_H) && ( (CyclinB && Cdk1) || (CyclinA && !Wee1 && Cdc25A) ) ;
    ge[Cdc25C_n] = f4N_DNA && Plk1 && ( (Cdc25B && !CHK1) || (CyclinB && Cdk1) ) ;
    ge[Cdk1_n] = CyclinB && Cdc25C && ( !CHK1 || (!Wee1 && Cdk1) ) ;
    ge[pAPC_n] = (CyclinB && Cdk1 && Plk1) || (CyclinB && Cdk1 && pAPC) || (pAPC && Cdc20) ;
    ge[Cdc20_n] = pAPC && !Emi1 && !Cdh1 && ( !Mad2 || ( !CyclinA && !(CyclinB && Cdk1) ) ) ;
    ge[Cdh1_n] = !(CyclinB && Cdk1) && !( CyclinA && (Emi1 || Cdc25A) ) ;
    ge[Replication_n] = !CAD && Pre_RC && ( (E2F1 && CyclinE && Cdc25A) || (Replication && CyclinA && Cdc25A && (E2F1 || !f4N_DNA) ))  ;
    ge[f4N_DNA_n] = !CAD && ( (Replication && ((Pre_RC && CyclinA) || f4N_DNA) ) || (f4N_DNA && !Ect2));
    ge[U_Kinetochores_n] = f4N_DNA && !Cdh1 && !A_Kinetochores && ( (CyclinB && Cdk1 && Cdc25C) || U_Kinetochores ) ;
    ge[Mad2_n] = U_Kinetochores && !A_Kinetochores ;
    ge[A_Kinetochores_n] = f4N_DNA && !Cdh1 && !(pAPC && Cdc20) && ( A_Kinetochores || (U_Kinetochores && Plk1 && CyclinB && Cdk1) ) ;
    ge[Plk1_H_n] = Plk1 && FoxM1 && (Plk1_H || FoxO3 || FoxO1) ;
    ge[Ect2_n] = f4N_DNA && Plk1_H && Cdh1 && !U_Kinetochores && !A_Kinetochores ;
    ge[Casp8_n] = DR4_5 || Casp3 ;
    ge[Casp2_n] = Casp3 || (U_Kinetochores && Mad2 && !(CyclinB && Cdk1) ) ;
    ge[MCL_1_n] = !Casp3 && !Casp2 && (!GSK3 || (AKT_B && (ERK || !E2F1))) && (!(Cdk1 && CyclinB && U_Kinetochores)) ;
    ge[BCLXL_n] = !Casp3 && (BCL2 || !BAD) && ( !U_Kinetochores || ( Plk1 && ( !(CyclinB && Cdk1) || (BCL2 && MCL_1) )) ||  ((BCL2 && MCL_1) && !(CyclinB && Cdk1) ) ) ;
    ge[BCL2_n] = !(Casp3 || BAD || BIM || BIK) &&  ( !U_Kinetochores || (MCL_1 && BCLXL) || (Plk1 && (BCLXL || MCL_1 || !(Cdk1 && CyclinB) ) ) ) ;
    ge[BAD_n] = Casp3 || !(AKT_H || AKT_B || ERK || S6K) || ( Casp8 && (!(AKT_B && ERK && S6K) && !(AKT_H && (AKT_B || ERK)) )) ;
    ge[BIK_n] = !(MCL_1 || BCLXL || BCL2) ;
    ge[BIM_n] = FoxO3 && GSK3 && !(ERK || MCL_1 || BCLXL || BCL2) ;
    ge[BID_n] = Casp8 || ( Casp2 && !(BCL2 || BCLXL || MCL_1) ) ;
    ge[BAK_n] = (BID && (BIM || BIK || !(BCL2 && BCLXL && MCL_1) ) ) || ((BIM || BIK) && !(BCLXL || MCL_1)) ;
    ge[BAX_n] = (BIM && ((BID || BIK) || !(BCL2 && BCLXL && MCL_1) ) ) || ((BID || BIK) && !(BCL2 || BCLXL)) ;
    ge[Cyto_C_n] = BAX || BAK ;
    ge[SMAC_n] = BAX || BAK ;
    ge[IAPs_n] = !SMAC || AKT_H ;
    ge[Casp9_n] = Casp3 || (!IAPs && Cyto_C) ;
    ge[Casp3_n] = (Casp9 && Casp8) || ( Casp3 && (Casp9 || Casp8) ) || ( !IAPs && (Casp9 || Casp8 || Casp3) ) ;
    ge[CAD_n] = Casp3 && Casp9 ;

}


unsigned int check_gates(Boolean_Dynamics_TRACKER *A){
    unsigned int *sin_now;
    int ok=1;
    snode::slinklist::slink *sc2;
    
    
    ge=new unsigned int [A->BRN->N+1];
    for (unsigned int i=1; i<=A->BRN->N+1; i++) {
        ge[i]=-1;
    }
    
    for (unsigned int ind1=1;ind1<=A->BRN->N;ind1++){
        ok=1;
        sc2=A->BRN->BN->node[ind1].Li.first_link_pt;
        while(sc2!=NULL){
            printf("%s\t\t",A->BRN->MDAT->Node_name[sc2->sneighbor]);
            sc2=sc2->next(sc2);
        }
        printf(">> %s\t\n",A->BRN->MDAT->Node_name[ind1]);

        sin_now=new unsigned int[A->BRN->N+1];
        for(unsigned long long k=0;k<A->BRN->Gate[ind1]->N;k++){
            decimal_to_any_base(k,2,A->BRN->Gate[ind1]->K,sin_now);
            
            sc2=A->BRN->BN->node[ind1].Li.first_link_pt;
            
            for(unsigned long long l=1;l<=A->BRN->Gate[ind1]->K;l++){
                printf("%d\t\t",sin_now[A->BRN->Gate[ind1]->K-l]);
                A->s[sc2->sneighbor]=sin_now[A->BRN->Gate[ind1]->K-l];
                sc2=sc2->next(sc2);
            }
            printf("\t\t | %d\n",A->BRN->Gate[ind1]->O[k]);
            //A->update_Gate(ind1); // should be the same as the gate at [k] !
            evaluate_gate_expressions(A);
            if (ge[ind1]==-1) {
                printf("No logic gate for node %s yet\n",A->BRN->MDAT->Node_name[ind1]);
            }
            else {
            if(ge[ind1] != A->BRN->Gate[ind1]->O[k]) {
                ok=0;
                printf("Gate %s row %lld does not match: table=%d expr=%d\n",A->BRN->MDAT->Node_name[ind1],k+1,A->BRN->Gate[ind1]->O[k],ge[ind1]);
                }
            }
        }
        if(ok==0) {
            printf("Logic gate for node %s is incorrect\n",A->BRN->MDAT->Node_name[ind1]);
            getchar();
        }
        delete[] sin_now;sin_now=NULL;
    }
    return(ok);
}

/*
unsigned short int FoxO3(unsigned short int AKT_B, unsigned short int AKT_H,unsigned short int ERK, unsigned short int Plk1){
    a = (!(AKT_B || AKT_H || ERK) || ( (!AKT_H || nor (ERK && AKT_B)) && !(Plk1 || Plk1_H));
    return(a);
}
  */


#endif /* Verify_gates_h */
