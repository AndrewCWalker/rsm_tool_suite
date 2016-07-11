#ifndef RSM_H
#define RSM_H

double CdSetup(struct initCd_struct *initCd, struct RSM_struct *RSMdata);

void emuInit(struct RSM_struct *RSMdata);

void emu(struct RSM_struct *RSMdata, double *xstar, double *specw, double *ystar);

double AreaLookUp(struct area_struct *AreaData, struct initCd_struct *initCd);

void read_input_file(struct initCd_struct *initCd);

void read_RSM(struct initCd_struct *initCd, struct RSM_struct *RSMdata);

void read_area(struct initCd_struct *initCd, struct area_struct *AreaData);

double BilinearInterp(double *pitch_array, double *yaw_array, double pitch, double yaw, double param[], int ipitch, int iyaw);

void WriteOutput(struct initCd_struct *initCd, double RSMCd, double RSMA, double RSMBC, char *adsstr);

#endif
