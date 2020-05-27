import configparser
import mysql.connector
import datetime
import math
import numpy as np
from decimal import Decimal
from pulp import *
import MySQLdb

# 24时刻太阳高度角
# 计算日在一年中的日期序号   输入时刻
def getRh(Ndate, Tm, cursor):
    # 时差修正
    cursor.execute(
        "select sc from taskinvoke_fhyc_scxz where id = %s ", (Ndate,))
    sc = float(cursor.fetchall()[0][0])
    delta = 23.45 * math.sin(360 / 180 * math.pi * (284 + Ndate) / 365) / 180 * math.pi  # 赤纬角计算公式1（弧度）
    Lcc = 115  # 保定地区当地子午线经度
    Rfai = 39 / 180 * math.pi  # 保定地区纬度，以弧度表示
    Lm = 120  # 东八区中央子午线经度
    RT = Tm + (Lcc - Lm) / 15 + sc / 60  # 当地的真太阳时
    Rbeta = (RT - 12) * 15 / 180 * math.pi  # 太阳时角
    Rh = math.asin(math.cos(Rfai) * math.cos(delta) * math.cos(Rbeta) + math.sin(Rfai) * math.sin(delta))  # 24时刻太阳高度角
    return Rh


# 水平面太阳总辐射照度 W/m2
def getIhdata(Ndate, Tm, date, cityid, cursor):
    Po = getRh(Ndate, Tm, cursor)
    P = 0
    if Po > 0:
        P = 1
    cursor.execute('select tqqk from taskinvoke_tq_hour where cityid = %s and ybrq = %s and hour = %s',
                   (cityid, date, Tm))
    tqqk = cursor.fetchall()[0][0]
    xs = 0
    if tqqk == '晴':
        xs = 1.12
    elif tqqk.find('雨') != -1 or tqqk.find('雪') != -1:
        xs = 0.44
    else:
        xs = 0.91

    zd = 'abc'
    if date.month == 11:
        zd = 'fsqd11'
    elif date.month == 12:
        if date.day < 15:
            zd = 'fsqd12s'
        elif date.day >= 15:
            zd = 'fsqd12x'
    elif date.month == 1:
        if date.day < 15:
            zd = 'fsqd1s'
        elif date.day >= 15:
            zd = 'fsqd1x'
    elif date.month == 2:
        if date.day < 15:
            zd = 'fsqd2s'
        elif date.day >= 15:
            zd = 'fsqd2x'
    elif date.month == 3:
        zd = 'fsqd3'
    sql = ''
    if zd == 'abc':
        sql = "select fssjd from taskinvoke_tq_tyfslbsj where fslb = 1 "
    else:
        sql = "select fssjd," + zd + " from taskinvoke_tq_tyfslbsj where fslb = 1 "
    cursor.execute(sql)
    a = cursor.fetchall()
    weather = 0
    for i in a:
        if len(i) == 2:
            if Tm == i[0] or Tm == int(i[0]):
                weather = float(i[1])
    return weather * xs * P


def getYx(lsxgid, nyzid, sblxid, cursor):
    cursor.execute(
        'select cysbid from taskinvoke_yhdd_clcysbls where lsxgid = %s and nyzid = %s and sblxid = %s',
        (lsxgid, nyzid, sblxid))
    dyrbryyx = cursor.fetchall()
    dyrbyx = 0  # 地源热泵预选
    if len(dyrbryyx) != 0:
        dyrbyx = 1
    return dyrbyx


def getTs(lsxgid, nyzid, sblxid, cursor):
    cursor.execute(
        'select cysbid from taskinvoke_yhdd_clcysbls where lsxgid = %s and nyzid = %s and sblxid = %s',
        (lsxgid, nyzid, sblxid))
    dyrb = cursor.fetchall()
    dyrbts = [0, 0, 0, 0, 0, 0]
    if len(dyrb) != 0:
        dyy = []
        for dy in range(len(dyrb)):
            dyy.append(1)
        if len(dyy) != 6:
            for dd in range(6 - len(dyy)):
                dyy.append(0)
        dyrbts = dyy
    return dyrbts


def getCs(cs1, cs2, cs3, cs4, lsxgid, nyzid, sblxid, cursor, nyzlx):
    sql = 'select ' + cs1 + ',' + cs2 + ',' + cs3 + ',' + cs4 + ' from taskinvoke_yhdd_clcysbls where lsxgid = %s and nyzid = %s and sblxid = %s'
    cursor.execute(sql, (lsxgid, nyzid, sblxid))
    dyrbcs = cursor.fetchall()
    dyrbrl = []
    dyrbxl = []
    dyrbzdfhl = []
    dyrbsbhdxs = []
    for cs in dyrbcs:
        if isinstance(cs[0], Decimal):
            dyrbrl.append(float(cs[0]))
        else:
            dyrbrl.append(cs[0])
        if isinstance(cs[1], Decimal):
            dyrbxl.append(float(cs[1]))
        else:
            dyrbxl.append(cs[1])
        if isinstance(cs[2], Decimal):
            dyrbzdfhl.append(float(cs[2]))
        else:
            dyrbzdfhl.append(cs[2])
        if isinstance(cs[3], Decimal):
            dyrbsbhdxs.append(float(cs[3]))
        else:
            dyrbsbhdxs.append(cs[3])
    if len(dyrbrl) != 6:
        for dd in range(6 - len(dyrbrl)):
            dyrbrl.append(0)
    if len(dyrbxl) != 6:
        for dd in range(6 - len(dyrbxl)):
            if sblxid == '01':
                dyrbxl.append(4.5)
            else:
                dyrbxl.append(0)
    if len(dyrbzdfhl) != 6:
        for dd in range(6 - len(dyrbzdfhl)):
            if sblxid == '01':
                dyrbzdfhl.append(0.3)
            else:
                dyrbzdfhl.append(0)
    if len(dyrbsbhdxs) != 6:
        for dd in range(6 - len(dyrbsbhdxs)):
            if sblxid == '01':
                dyrbsbhdxs.append(0.2)
            elif sblxid == '02':
                dyrbsbhdxs.append(0.3)
    list = []
    list.append(dyrbrl)
    list.append(dyrbxl)
    list.append(dyrbzdfhl)
    list.append(dyrbsbhdxs)
    return list


def getCs2(cs1, cs2, cs3, lsxgid, nyzid, sblxid, cursor, nyzlx):
    sql = 'select ' + cs1 + ',' + cs2 + ',' + cs3 + ' from taskinvoke_yhdd_clcysbls where lsxgid = %s and nyzid = %s and sblxid = %s'
    cursor.execute(sql, (lsxgid, nyzid, sblxid))
    dyrbcs = cursor.fetchall()
    dyrbrl = []
    dyrbxl = []
    dyrbzdfhl = []
    for cs in dyrbcs:
        if isinstance(cs[0], Decimal):
            dyrbrl.append(float(cs[0]))
        else:
            dyrbrl.append(cs[0])
        if isinstance(cs[1], Decimal):
            dyrbxl.append(float(cs[1]))
        else:
            dyrbxl.append(cs[1])
        if isinstance(cs[2], Decimal):
            dyrbzdfhl.append(float(cs[2]))
        else:
            dyrbzdfhl.append(cs[2])
    if len(dyrbrl) != 6:
        for dd in range(6 - len(dyrbrl)):
            dyrbrl.append(0)
    if len(dyrbxl) != 6:
        for dd in range(6 - len(dyrbxl)):
            if nyzlx == '1':
                dyrbxl.append(0)
            else:
                if sblxid == '03':
                    dyrbxl.append(0.89)
                elif sblxid == '04':
                    dyrbxl.append(0.98)
    if len(dyrbzdfhl) != 6:
        for dd in range(6 - len(dyrbzdfhl)):
            if sblxid == '03':
                dyrbzdfhl.append(0.26)
            elif sblxid == '04':
                dyrbzdfhl.append(0.1)
    list = []
    list.append(dyrbrl)
    list.append(dyrbxl)
    list.append(dyrbzdfhl)
    return list


def get(preselection_regi, preselection_GSHP, preselection_GICE, preselection_GB, preselection_EB,
        GSHP_select_configured, GSHP_capacity, GSHP_minplr, eff_GSHP, GSHP_pump_ele,
        GICE_select_configured, GICE_rated_power, effE_GICE, effH_GICE, GICE_minplr,
        GB_select_configured, GB_capacity, eff_GB, GB_minplr,
        EB_select_configured, EB_capacity, eff_EB, EB_minplr,
        preselection_block, preselection_GSHP_X, preselection_GICE_X, preselection_GB_X, preselection_EB_X, preselection_HC_X, preselection_WTSH_X,
        GSHP_select_configured_X, GSHP_capacity_X, eff_GSHP_X, GSHP_minplr_X, GSHP_pump_ele_X,
        GICE_select_configured_X, GICE_rated_power_X, effE_GICE_X, effH_GICE_X, GICE_minplr_X,
        GB_select_configured_X, GB_capacity_X, eff_GB_X, GB_minplr_X,
        EB_select_configured_X, EB_capacity_X, eff_EB_X, EB_minplr_X,
        HC_area_X, eff_HC_X,
        WTSH_voulme_X, eff_Hin_X, eff_Hout_X, WTSH_loss_X, ratio_Hin_X, ratio_Hout_X,
        eff_pipenet_X, eff_HE_X,
        preselection_user, preselection_GSHP_Y, preselection_GICE_Y, preselection_GB_Y, preselection_EB_Y, preselection_HC_Y, preselection_WTSH_Y,
        GSHP_select_configured_Y, GSHP_capacity_Y, eff_GSHP_Y, GSHP_minplr_Y, GSHP_pump_ele_Y,
        GICE_select_configured_Y, GICE_rated_power_Y, effE_GICE_Y, effH_GICE_Y, GICE_minplr_Y,
        GB_select_configured_Y, GB_capacity_Y, eff_GB_Y, GB_minplr_Y,
        EB_select_configured_Y, EB_capacity_Y, eff_EB_Y, EB_minplr_Y,
        HC_area_Y, eff_HC_Y,
        WTSH_voulme_Y, eff_Hin_Y, eff_Hout_Y, WTSH_loss_Y, ratio_Hin_Y, ratio_Hout_Y, eff_pipenet_Y,
        Load_hea, I_solar, Ele_price, Ele_sell_price, Gas_price):
    hh = list(range(24))
    xx = list(range(13))
    yy = list(range(3))
    ii = list(range(6))
    # 区域能源站
    # 地源热泵小时出力、小时耗电量、启停变量
    GSHP_output = LpVariable.dicts('GSHP_output', ((i, h) for i in ii for h in hh), lowBound=0)
    GSHP_ele_comumption = LpVariable.dicts('GSHP_ele_comumption', ((i, h) for i in ii for h in hh), lowBound=0)
    GSHP_binary = LpVariable.dicts('GSHP_binary', ((i, h) for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 燃气内燃机小时发电量、余热回收量、燃气消耗量、售电量、启停变量
    GICE_ele_output = LpVariable.dicts('GICE_ele_output', ((i, h) for i in ii for h in hh), lowBound=0)
    GICE_hea_output = LpVariable.dicts('GICE_hea_output', ((i, h) for i in ii for h in hh), lowBound=0)
    GICE_gas_comsuption = LpVariable.dicts('GICE_gas_comsuption', ((i, h) for i in ii for h in hh))
    GICEEle_sell_grid = LpVariable.dicts('GICEEle_sell_grid', (h for h in hh), lowBound=0)
    GICE_binary = LpVariable.dicts('GICE_binary', ((i, h) for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 燃气锅炉小时出力、燃气消耗、启停变量
    GB_hea_output = LpVariable.dicts('GB_hea_output', ((i, h) for i in ii for h in hh), lowBound=0)
    GB_gas_comsuption = LpVariable.dicts('GB_gas_comsuption', ((i, h) for i in ii for h in hh), lowBound=0)
    GB_binary = LpVariable.dicts('GB_binary', ((i, h) for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 电锅炉小时出力、耗电量、启停变量
    EB_hea_output = LpVariable.dicts('EB_hea_output', ((i, h) for i in ii for h in hh), lowBound=0)
    EB_ele_comsuption = LpVariable.dicts('EB_ele_comsuption', ((i, h) for i in ii for h in hh), lowBound=0)
    EB_binary = LpVariable.dicts('EB_binary', ((i, h) for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)

    # 街区能源站
    # 地源热泵小时出力、小时耗电量、启停变量
    GSHP_output_X = LpVariable.dicts('GSHP_output_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0)
    GSHP_ele_comumption_X = LpVariable.dicts('GSHP_ele_comumption_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0)
    GSHP_binary_X = LpVariable.dicts('GSHP_binary_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 燃气内燃机小时发电量、余热回收量、燃气消耗量、售电量、启停变量
    GICE_ele_output_X = LpVariable.dicts('GICE_ele_output_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0)
    GICE_hea_output_X = LpVariable.dicts('GICE_hea_output_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0)
    GICE_gas_comsuption_X = LpVariable.dicts('GICE_gas_comsuption_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0)
    GICEEle_sell_grid_X = LpVariable.dicts('GICEEle_sell_grid_X', ((x, h) for x in xx for h in hh), lowBound=0)
    GICE_binary_X = LpVariable.dicts('GICE_binary_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 燃气锅炉小时出力、燃气消耗、启停变量
    GB_hea_output_X = LpVariable.dicts('GB_hea_output_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0)
    GB_gas_comsuption_X = LpVariable.dicts('GB_gas_comsuption_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0)
    GB_binary_X = LpVariable.dicts('GB_binary_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 电锅炉小时出力、耗电量、启停变量
    EB_hea_output_X = LpVariable.dicts('EB_hea_output_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0)
    EB_ele_comsuption_X = LpVariable.dicts('EB_ele_comsuption_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0)
    EB_binary_X = LpVariable.dicts('EB_binary_X', ((x, i, h) for x in xx for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 光热锅炉小时出力
    HC_output_X = LpVariable.dicts('HC_output_X', ((x, h) for x in xx for h in hh), lowBound=0)
    # 水蓄能小时储能、释能、储能量
    WTSH_in_X = LpVariable.dicts('WTSH_in_X', ((x, h) for x in xx for h in hh), lowBound=0)
    WTSH_out_X = LpVariable.dicts('WTSH_out_X', ((x, h) for x in xx for h in hh), lowBound=0)
    WTSH_X = LpVariable.dicts('WTSH_X', ((x, h) for x in xx for h in hh), lowBound=0)
    # 热交换器
    HE_hea_output_X = LpVariable.dicts('HE_hea_output_X', ((x, h) for x in xx for h in hh), lowBound=0)
    HE_binary_X = LpVariable.dicts('HE_binary_X', ((x, h) for x in xx for h in hh), lowBound=0)

    # 用户能源站
    # 地源热泵小时出力、小时耗电量、启停变量
    GSHP_output_Y = LpVariable.dicts('GSHP_output_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0)
    GSHP_ele_comumption_Y = LpVariable.dicts('GSHP_ele_comumption_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0)
    GSHP_binary_Y = LpVariable.dicts('GSHP_binary_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 燃气内燃机小时发电量、余热回收量，燃气消耗量、售电量、启停变量
    GICE_ele_output_Y = LpVariable.dicts('GICE_ele_output_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0)
    GICE_hea_output_Y = LpVariable.dicts('GICE_hea_output_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0)
    GICE_gas_comsuption_Y = LpVariable.dicts('GICE_gas_comsuption_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0)
    GICEEle_sell_grid_Y = LpVariable.dicts('GICEEle_sell_grid_Y', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)
    GICE_binary_Y = LpVariable.dicts('GICE_binary_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 燃气锅炉小时出力、燃气消耗、启停变量
    GB_hea_output_Y = LpVariable.dicts('GB_hea_output_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0)
    GB_gas_comsuption_Y = LpVariable.dicts('GB_gas_comsuption_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0)
    GB_binary_Y = LpVariable.dicts('GB_binary_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 电锅炉小时出力、耗电量
    EB_hea_output_Y = LpVariable.dicts('EB_hea_output_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0)
    EB_ele_comsuption_Y = LpVariable.dicts('EB_ele_comsuption_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0)
    EB_binary_Y = LpVariable.dicts('EB_binary_Y', ((x, y, i, h) for x in xx for y in yy for i in ii for h in hh), lowBound=0,upBound=1,cat=LpBinary)
    # 光热锅炉小时出力
    HC_output_Y = LpVariable.dicts('HC_output_Y', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)
    # 水蓄能小时储能、释能、储能量
    WTSH_in_Y = LpVariable.dicts('WTSH_in_Y', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)
    WTSH_out_Y = LpVariable.dicts('WTSH_out_Y', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)
    WTSH_Y = LpVariable.dicts('WTSH_Y', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)

    # 总运行费用、总购电费用、总燃气费用、总售电费用
    Tatal_cost = LpVariable('Tatal_cost')
    Ele_cost = LpVariable('Ele_cost')
    Gas_cost = LpVariable('Gas_cost')
    Ele_sell = LpVariable('Ele_sell')
    # 小时购电量、小时购燃气量、小时售电量
    Ele_totalbuy_grid = LpVariable.dicts('Ele_totalbuy_grid', (h for h in hh), lowBound=0)
    Gas_total_comumption = LpVariable.dicts('Gas_total_comumption', (h for h in hh), lowBound=0)
    GICETotalEle_sell_grid = LpVariable.dicts('GICETotalEle_sell_grid', (h for h in hh), lowBound=0)
    # 区域能源站购电量、街区能源站购电量、用户能源站购电量 h  x,h  x,y,h 123
    Ele_buy_grid = LpVariable.dicts('Ele_buy_grid', (h for h in hh), lowBound=0)
    Ele_buy_grid_X = LpVariable.dicts('Ele_buy_grid_X', ((x, h) for x in xx for h in hh), lowBound=0)
    Ele_buy_grid_Y = LpVariable.dicts('Ele_buy_grid_Y', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)
    # 区域能源站输给第x个街区能源站的热量、经管网输给第x个街区能源站的热量、区域能源站输出的总热量1
    input_block = LpVariable.dicts('input_block', ((x, h) for x in xx for h in hh), lowBound=0)
    input_block_X = LpVariable.dicts('input_block_X', ((x, h) for x in xx for h in hh), lowBound=0)
    output_regi = LpVariable.dicts('output_regi', (h for h in hh), lowBound=0)
    # 第x个街区能源站输给第y个用户能源站的热量、经管网输给第y个用户能源站的热量、第x个街区能源站输出的总热量1
    input_user = LpVariable.dicts('input_user', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)
    input_user_Y = LpVariable.dicts('input_user_Y', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)
    output_block = LpVariable.dicts('output_block', ((x, h) for x in xx for h in hh), lowBound=0)
    # 区域能源站输送给第x个街区能源站的电量, 区域能源站输出的总电量1
    input_ele_block = LpVariable.dicts('input_ele_block', ((x, h) for x in xx for h in hh), lowBound=0)
    output_ele_regi = LpVariable.dicts('output_ele_regi', (h for h in hh), lowBound=0)
    # 街区能源站输送给第y个街区能源站的电量、街区能源站输出的总电量1
    input_ele_user = LpVariable.dicts('input_ele_user', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)
    output_ele_block = LpVariable.dicts('output_ele_block', ((x, h) for x in xx for h in hh), lowBound=0)

    # 总购电费用  (每小时的购电价*每小时的购电量 之和)
    # 每小时从电网购电总量=区域能源站购电量+街区能源站购电量+用户能源站购电量：
    for h in hh:
        Ele_totalbuy_grid[h] = Ele_buy_grid[h] + sum(Ele_buy_grid_X[(x,h)] for x in xx) + sum(Ele_buy_grid_Y[(x,y,h)] for x in xx for y in yy)
    Ele_cost = sum(Ele_price[h] * Ele_totalbuy_grid[h] for h in hh)
    # 每小时燃气总消耗量=区域能源站+街区能源站+用户能源站
    for h in hh:
        Gas_total_comumption[h] = sum(GICE_gas_comsuption[(i,h)] + GB_gas_comsuption[(i, h)] for i in ii) + sum(GICE_gas_comsuption_X[(x, i, h)] + GB_gas_comsuption_X[(x, i, h)] for x in xx for i in ii) + sum(GICE_gas_comsuption_Y[(x, y, i, h)] + GB_gas_comsuption_Y[(x, y, i, h)] for x in xx for y in yy for i in ii)
    # 总购燃气费用  (燃气价格*24小时的燃气总消耗)
    Gas_cost = sum(Gas_total_comumption[h] for h in hh) * Gas_price
    # 每小时售电量=区域能源站售电量+街区能源站售电量+用户能源站售电量：
    for h in hh:
        GICETotalEle_sell_grid[h] = GICEEle_sell_grid[h] + sum(GICEEle_sell_grid_X[(x,h)] for x in xx) + sum(GICEEle_sell_grid_Y[(x, y,h)] for x in xx for y in yy)
    # 总售电费用  (售电价格*24小时的售电量)
    Ele_sell = sum(GICETotalEle_sell_grid[h] for h in hh)* Ele_sell_price
    # 日总运行费用
    Tatal_cost = Ele_cost + Gas_cost - Ele_sell

    # 创建问题实例，求最小极值
    prob = LpProblem("YHDD_riqian_Problem", LpMinimize)
    # 添加目标方程
    prob += Tatal_cost

    # 添加约束条件
    # 区域能源站电力实时平衡约束
    # 电网购电量 + 燃气内燃机发电量 >= 地源热泵耗电量 + 电锅炉耗电量 + 燃气内燃机售电量 + 向街区能源站输送的总电量
    for h in range(24):
        prob += Ele_buy_grid[h] + GICE_ele_output[(0,h)] + GICE_ele_output[(1,h)] + GICE_ele_output[(2,h)] + GICE_ele_output[(3,h)] + GICE_ele_output[(4,h)] + GICE_ele_output[(5,h)] >= GSHP_ele_comumption[(0,h)] + GSHP_ele_comumption[(1,h)] + GSHP_ele_comumption[(2,h)] + GSHP_ele_comumption[(3,h)] + GSHP_ele_comumption[(4,h)] + GSHP_ele_comumption[(5,h)] + EB_ele_comsuption[(0,h)] + EB_ele_comsuption[(1,h)] + EB_ele_comsuption[(2,h)] + EB_ele_comsuption[(3,h)] + EB_ele_comsuption[(4,h)] + EB_ele_comsuption[(5,h)] + GICEEle_sell_grid[h] + input_ele_block[(0,h)] + input_ele_block[(1,h)] + input_ele_block[(2,h)] + input_ele_block[(3,h)] + input_ele_block[(4,h)] + input_ele_block[(5,h)] + input_ele_block[(6,h)] + input_ele_block[(7,h)] + input_ele_block[(8,h)] + input_ele_block[(9,h)] + input_ele_block[(10,h)] + input_ele_block[(11,h)] + input_ele_block[(12,h)]
        prob += GICE_ele_output[(0,h)] + GICE_ele_output[(1,h)] + GICE_ele_output[(2,h)] + GICE_ele_output[(3,h)] + GICE_ele_output[(4,h)] + GICE_ele_output[(5,h)] >= GICEEle_sell_grid[h] + input_ele_block[(0,h)] + input_ele_block[(1,h)] + input_ele_block[(2,h)] + input_ele_block[(3,h)] + input_ele_block[(4,h)] + input_ele_block[(5,h)] + input_ele_block[(6,h)] + input_ele_block[(7,h)] + input_ele_block[(8,h)] + input_ele_block[(9,h)] + input_ele_block[(10,h)] + input_ele_block[(11,h)] + input_ele_block[(12,h)]
    # 区域能源站热力实时平衡
    # 地源热泵产热量 + 燃气内燃机产热量 + 燃气锅炉产热 + 电锅炉产热 >= 送给各个街区能源站的热量总和
    for h in range(24):
        prob += GSHP_output[(0,h)] + GSHP_output[(1,h)] + GSHP_output[(2,h)] + GSHP_output[(3,h)] + GSHP_output[(4,h)] + GSHP_output[(5,h)] + GICE_hea_output[(0,h)] + GICE_hea_output[(1,h)] + GICE_hea_output[(2,h)] + GICE_hea_output[(3,h)] + GICE_hea_output[(4,h)] + GICE_hea_output[(5,h)] + GB_hea_output[(0,h)] + GB_hea_output[(1,h)] + GB_hea_output[(2,h)] + GB_hea_output[(3,h)] + GB_hea_output[(4,h)] + GB_hea_output[(5,h)] + EB_hea_output[(0,h)] + EB_hea_output[(1,h)] + EB_hea_output[(2,h)] + EB_hea_output[(3,h)] + EB_hea_output[(4,h)] + EB_hea_output[(5,h)] >= input_block[(0,h)] + input_block[(1,h)] + input_block[(2,h)] + input_block[(3,h)] + input_block[(4,h)] + input_block[(5,h)] + input_block[(6,h)] + input_block[(7,h)] + input_block[(8,h)] + input_block[(9,h)] + input_block[(10,h)] + input_block[(11,h)] + input_block[(12,h)]
    # 第x个街区能源站电力实时平衡
    # 电网购电量 + 燃气内燃机发电量 + 区域能源站输入电量 >= 地源热泵耗电量 + 电锅炉耗电量 + 燃气内燃机售电量 + 向用户能源站输送的总电量
    for x in range(13):
        for h in range(24):
            prob += Ele_buy_grid_X[(x,h)] + GICE_ele_output_X[(x,0,h)] + GICE_ele_output_X[(x,1,h)] + GICE_ele_output_X[(x,2,h)] + GICE_ele_output_X[(x,3,h)] + GICE_ele_output_X[(x,4,h)] + GICE_ele_output_X[(x,5,h)] + input_ele_block[(x,h)] >= GSHP_ele_comumption_X[(x,0,h)] + GSHP_ele_comumption_X[(x,1,h)] + GSHP_ele_comumption_X[(x,2,h)] + GSHP_ele_comumption_X[(x,3,h)] + GSHP_ele_comumption_X[(x,4,h)] + GSHP_ele_comumption_X[(x,5,h)] + EB_ele_comsuption_X[(x,0,h)] + EB_ele_comsuption_X[(x,1,h)] + EB_ele_comsuption_X[(x,2,h)] + EB_ele_comsuption_X[(x,3,h)] + EB_ele_comsuption_X[(x,4,h)] + EB_ele_comsuption_X[(x,5,h)] + GICEEle_sell_grid_X[(x,h)] + input_ele_user[(x,0,h)] + input_ele_user[(x,1,h)] + input_ele_user[(x,2,h)]
            prob += GICE_ele_output_X[(x,0,h)] + GICE_ele_output_X[(x,1,h)] + GICE_ele_output_X[(x,2,h)] + GICE_ele_output_X[(x,3,h)] + GICE_ele_output_X[(x,4,h)] + GICE_ele_output_X[(x,5,h)] + input_ele_block[(x,h)] >= GICEEle_sell_grid_X[(x,h)] + input_ele_user[(x,0,h)] + input_ele_user[(x,1,h)] + input_ele_user[(x,2,h)]
    # 第x个街区能源站热力实时平衡
    # 地源热泵产热量 + 燃气内燃机产热量 + 燃气锅炉产热 + 电锅炉产热 + 光热 + 释热量 - 蓄热量 + 区域能源站输入热量 + 街区能源站输入的电量 >= 输入用户能源站热量之和
    for x in range(13):
        for h in range(24):
            prob += GSHP_output_X[(x,0,h)] + GSHP_output_X[(x,1,h)] + GSHP_output_X[(x,2,h)] + GSHP_output_X[(x,3,h)] + GSHP_output_X[(x,4,h)] + GSHP_output_X[(x,5,h)] + GICE_hea_output_X[(x,0,h)] + GICE_hea_output_X[(x,1,h)] + GICE_hea_output_X[(x,2,h)] + GICE_hea_output_X[(x,3,h)] + GICE_hea_output_X[(x,4,h)] + GICE_hea_output_X[(x,5,h)] + GB_hea_output_X[(x,0,h)] + GB_hea_output_X[(x,1,h)] + GB_hea_output_X[(x,2,h)] + GB_hea_output_X[(x,3,h)] + GB_hea_output_X[(x,4,h)] + GB_hea_output_X[(x,5,h)] + EB_hea_output_X[(x,0,h)] + EB_hea_output_X[(x,1,h)] + EB_hea_output_X[(x,2,h)] + EB_hea_output_X[(x,3,h)] + EB_hea_output_X[(x,4,h)] + EB_hea_output_X[(x,5,h)] + HC_output_X[(x,h)] - WTSH_in_X[(x,h)] + WTSH_out_X[(x,h)] + HE_hea_output_X[(x,h)] >= input_user[(x,0,h)] + input_user[(x,1,h)] + input_user[(x,2,h)]
    # 第x个街区能源站第y个用户能源站电力实时平衡
    # 电网购电量 + 燃气内燃机发电量 >= 地源热泵耗电量 + 电锅炉耗电量 + 内燃机售电量
    for x in range(13):
        for y in range(3):
            for h in range(24):
                prob += Ele_buy_grid_Y[(x,y,h)] + GICE_ele_output_Y[(x,y,0,h)] + GICE_ele_output_Y[(x,y,1,h)] + GICE_ele_output_Y[(x,y,2,h)] + GICE_ele_output_Y[(x,y,3,h)] + GICE_ele_output_Y[(x,y,4,h)] + GICE_ele_output_Y[(x,y,5,h)] + input_ele_user[(x,y,h)] >= GSHP_ele_comumption_Y[(x,y,0,h)] + GSHP_ele_comumption_Y[(x,y,1,h)] + GSHP_ele_comumption_Y[(x,y,2,h)] + GSHP_ele_comumption_Y[(x,y,3,h)] + GSHP_ele_comumption_Y[(x,y,4,h)] + GSHP_ele_comumption_Y[(x,y,5,h)] + EB_ele_comsuption_Y[(x,y,0,h)] + EB_ele_comsuption_Y[(x,y,1,h)] + EB_ele_comsuption_Y[(x,y,2,h)] + EB_ele_comsuption_Y[(x,y,3,h)] + EB_ele_comsuption_Y[(x,y,4,h)] + EB_ele_comsuption_Y[(x,y,5,h)] +  GICEEle_sell_grid_Y[(x,y,h)]
    # 第x个街区能源站第y个用户能源站热力实时平衡
    # 地源热泵产热量 + 燃气内燃机产热量 + 燃气锅炉产热 + 电锅炉产热 + 光热 + 释热量 - 蓄热量 + 街区能源站输入热量 >= 热负荷
    for x in range(13):
        for y in range(3):
            for h in range(24):
                prob += GSHP_output_Y[(x,y,0,h)] + GSHP_output_Y[(x,y,1,h)] + GSHP_output_Y[(x,y,2,h)] + GSHP_output_Y[(x,y,3,h)] + GSHP_output_Y[(x,y,4,h)] + GSHP_output_Y[(x,y,5,h)] + GICE_hea_output_Y[(x,y,0,h)] + GICE_hea_output_Y[(x,y,1,h)] + GICE_hea_output_Y[(x,y,2,h)] + GICE_hea_output_Y[(x,y,3,h)] + GICE_hea_output_Y[(x,y,4,h)] + GICE_hea_output_Y[(x,y,5,h)] + GB_hea_output_Y[(x,y,0,h)] + GB_hea_output_Y[(x,y,1,h)] + GB_hea_output_Y[(x,y,2,h)] + GB_hea_output_Y[(x,y,3,h)] + GB_hea_output_Y[(x,y,4,h)] + GB_hea_output_Y[(x,y,5,h)] + EB_hea_output_Y[(x,y,0,h)] + EB_hea_output_Y[(x,y,1,h)] + EB_hea_output_Y[(x,y,2,h)] + EB_hea_output_Y[(x,y,3,h)] + EB_hea_output_Y[(x,y,4,h)] + EB_hea_output_Y[(x,y,5,h)] + HC_output_Y[(x,y,h)] - WTSH_in_Y[(x,y,h)] + WTSH_out_Y[(x,y,h)] + input_user[(x,y,h)] * eff_pipenet_Y[x][y] >= Load_hea[x][y][h]

    # 区域能源站
    # 区域能源站预选约束
    q_gas = 35880
    prob += preselection_GSHP >= 0
    prob += preselection_GSHP <= preselection_regi
    prob += preselection_GICE >= 0
    prob += preselection_GICE <= preselection_regi
    prob += preselection_GB >= 0
    prob += preselection_GB <= preselection_regi
    prob += preselection_EB >= 0
    prob += preselection_EB <= preselection_regi
    for i in range(6):
        # 地源热泵
        prob += GSHP_select_configured[i] >= 0
        prob += GSHP_select_configured[i] <= preselection_GSHP
        # 燃气内燃机
        prob += GICE_select_configured[i] >= 0
        prob += GICE_select_configured[i] <= preselection_GICE
        # 燃气锅炉
        prob += GB_select_configured[i] >= 0
        prob += GB_select_configured[i] <= preselection_GB
        # 电锅炉
        prob += EB_select_configured[i] >= 0
        prob += EB_select_configured[i] <= preselection_EB
        for h in range(24):
            # 地源热泵
            prob += GSHP_binary[(i,h)] * GSHP_capacity[i] * GSHP_minplr[i] <= GSHP_output[(i,h)]
            prob += GSHP_output[(i,h)] <= GSHP_capacity[i] * GSHP_binary[(i,h)]
            prob += GSHP_ele_comumption[(i,h)] == GSHP_output[(i,h)] * (1/ eff_GSHP[i]) * (1 + GSHP_pump_ele[i])
            prob += (GSHP_binary[(i,h + 2)] - GSHP_binary[(i,h + 1)]) - (GSHP_binary[(i,h + 1)] - GSHP_binary[(i,h)]) >= -1
            prob += (GSHP_binary[(i,h + 2)] - GSHP_binary[(i,h + 1)]) - (GSHP_binary[(i,h + 1)] - GSHP_binary[(i,h)]) <= 1
            # 燃气内燃机
            prob += GICE_binary[(i,h)] * GICE_rated_power[i] * GICE_minplr[i] <= GICE_ele_output[(i,h)]
            prob += GICE_ele_output[(i,h)] <= GICE_rated_power[i] * GICE_binary[(i,h)]
            prob += GICE_ele_output[(i,h)] == GICE_gas_comsuption[(i,h)] * effE_GICE[i] * q_gas / 3600
            prob += GICE_hea_output[(i,h)] == GICE_gas_comsuption[(i,h)] * effH_GICE[i] * q_gas / 3600
            prob += (GICE_binary[(i,h + 2)] - GICE_binary[(i,h + 1)]) - (GICE_binary[(i,h + 1)] - GICE_binary[(i,h)]) >= -1
            prob += (GICE_binary[(i,h + 2)] - GICE_binary[(i,h + 1)]) - (GICE_binary[(i,h + 1)] - GICE_binary[(i,h)]) <= 1
            # 燃气锅炉
            prob += GB_binary[(i,h)] * GB_capacity[i] * GB_minplr[i] <= GB_hea_output[(i,h)]
            prob += GB_hea_output[(i,h)] <= GB_capacity[i] * GB_binary[(i,h)]
            prob += GB_hea_output[(i,h)] == GB_gas_comsuption[(i,h)] * q_gas * eff_GB[i] / 3600
            prob += (GB_binary[(i,h + 2)] - GB_binary[(i,h + 1)]) - (GB_binary[(i,h + 1)] - GB_binary[(i,h)]) >= -1
            prob += (GB_binary[(i,h + 2)] - GB_binary[(i,h + 1)]) - (GB_binary[(i,h + 1)] - GB_binary[(i,h)]) <= 1
            # 电锅炉
            prob += EB_binary[(i,h)] * EB_capacity[i] * EB_minplr[i] <= EB_hea_output[(i,h)]
            prob += EB_hea_output[(i,h)] <= EB_capacity[i] * EB_binary[(i,h)]
            prob += EB_hea_output[(i,h)] == EB_ele_comsuption[(i,h)] * eff_EB[i]
            prob += (EB_binary[(i,h + 2)] - EB_binary[(i,h + 1)]) - (EB_binary[(i,h + 1)] - EB_binary[(i,h)]) >= -1
            prob += (EB_binary[(i,h + 2)] - EB_binary[(i,h + 1)]) - (EB_binary[(i,h + 1)] - EB_binary[(i,h)]) <= 1
    for h in range(24):
        # 地源热泵
        prob += GSHP_binary[(0,h)] + GSHP_binary[(1,h)] + GSHP_binary[(2,h)] + GSHP_binary[(3,h)] + GSHP_binary[(4,h)] + GSHP_binary[(5,h)] <= sum(GSHP_select_configured)
        # 燃气内燃机
        prob += GICE_binary[(0,h)] + GICE_binary[(1,h)] + GICE_binary[(2,h)] + GICE_binary[(3,h)] + GICE_binary[(4,h)] + GICE_binary[(5,h)] <= sum(GICE_select_configured)
        prob += GICE_ele_output[(0,h)] + GICE_ele_output[(1,h)] + GICE_ele_output[(2,h)] + GICE_ele_output[(3,h)] + GICE_ele_output[(4,h)] + GICE_ele_output[(5,h)] >= GICEEle_sell_grid[h]
        # 燃气锅炉
        prob += GB_binary[(0,h)] + GB_binary[(1,h)] + GB_binary[(2,h)] + GB_binary[(3,h)] + GB_binary[(4,h)] + GB_binary[(5,h)] <= sum(GB_select_configured)
        # 电锅炉
        prob += EB_binary[(0,h)] + EB_binary[(1,h)] + EB_binary[(2,h)] + EB_binary[(3,h)] + EB_binary[(4,h)] + EB_binary[(5,h)] <= sum(EB_select_configured)
    # 地源热泵
    prob += sum(GSHP_select_configured) <= 6
    # 燃气内燃机
    prob += sum(GICE_select_configured) <= 6
    # 燃气锅炉
    prob += sum(GB_select_configured) <= 6
    # 电锅炉
    prob += sum(EB_select_configured) <= 6

    # 街区能源站
    M = 9000000000000000000000000
    density_H = 0.98
    for x in range(13):
        # 街区能源站预选约束
        prob += preselection_GSHP_X[x] >= 0
        prob += preselection_GSHP_X[x] <= preselection_block[x]
        prob += preselection_GICE_X[x] >= 0
        prob += preselection_GICE_X[x] <= preselection_block[x]
        prob += preselection_GB_X[x] >= 0
        prob += preselection_GB_X[x] <= preselection_block[x]
        prob += preselection_EB_X[x] >= 0
        prob += preselection_EB_X[x] <= preselection_block[x]
        prob += preselection_HC_X[x] >= 0
        prob += preselection_HC_X[x] <= preselection_block[x]
        prob += preselection_WTSH_X[x] >= 0
        prob += preselection_WTSH_X[x] <= preselection_block[x]
        for i in range(6):
            # 地源热泵
            prob += GSHP_select_configured_X[x][i] >= 0
            prob += GSHP_select_configured_X[x][i] <= preselection_GSHP_X[x]
            # 燃气内燃机
            prob += GICE_rated_power_X[x][i] >= 0
            prob += GICE_rated_power_X[x][i] <= preselection_GICE_X[x]
            # 燃气锅炉
            prob += GB_select_configured_X[x][i] >= 0
            prob += GB_select_configured_X[x][i] <= preselection_GB_X[x]
            # 电锅炉
            prob += EB_select_configured_X[x][i] >= 0
            prob += EB_select_configured_X[x][i] <= preselection_EB_X[x]
            for h in range(24):
                # 地源热泵
                prob += GSHP_binary_X[(x,i,h)] * GSHP_capacity_X[x][i] * GSHP_minplr_X[x][i] <= GSHP_output_X[(x,i,h)]
                prob += GSHP_output_X[(x,i,h)] <= GSHP_capacity_X[x][i] * GSHP_binary_X[(x,i,h)]
                prob += GSHP_ele_comumption_X[(x,i,h)] == GSHP_output_X[(x,i,h)] * (1/ eff_GSHP_X[x][i]) * (1 + GSHP_pump_ele_X[x][i])
                prob += (GSHP_binary_X[(x,i,h + 2)] - GSHP_binary_X[(x,i,h + 1)]) - (GSHP_binary_X[(x,i,h + 1)] - GSHP_binary_X[(x,i,h)]) >= -1
                prob += (GSHP_binary_X[(x,i,h + 2)] - GSHP_binary_X[(x,i,h + 1)]) - (GSHP_binary_X[(x,i,h + 1)] - GSHP_binary_X[(x,i,h)]) <= 1
                # 燃气内燃机
                prob += GICE_binary_X[(x,i,h)] * GICE_rated_power_X[x][i] * GICE_minplr_X[x][i] <= GICE_ele_output_X[(x,i,h)]
                prob += GICE_ele_output_X[(x,i,h)] <= GICE_rated_power_X[x][i] * GICE_binary_X[(x,i,h)]
                prob += GICE_ele_output_X[(x,i,h)] == GICE_gas_comsuption_X[(x,i,h)] * effE_GICE_X[x][i] * q_gas / 3600
                prob += GICE_hea_output_X[(x,i,h)] == GICE_gas_comsuption_X[(x,i,h)] * effH_GICE_X[x][i] * q_gas / 3600
                prob += (GICE_binary_X[(x,i,h + 2)] - GICE_binary_X[(x,i,h + 1)]) - (GICE_binary_X[(x,i,h + 1)] - GICE_binary_X[(x,i,h)]) >= -1
                prob += (GICE_binary_X[(x,i,h + 2)] - GICE_binary_X[(x,i,h + 1)]) - (GICE_binary_X[(x,i,h + 1)] - GICE_binary_X[(x,i,h)]) <= 1
                # 燃气锅炉
                prob += GB_binary_X[(x,i,h)] * GB_capacity_X[x][i] * GB_minplr_X[x][i] <= GB_hea_output_X[(x,i,h)]
                prob += GB_hea_output_X[(x,i,h)] <= GB_capacity_X[x][i] * GB_binary_X[(x,i,h)]
                prob += GB_hea_output_X[(x,i,h)] == GB_gas_comsuption_X[(x,i,h)] * q_gas * eff_GB_X[x][i] / 3600
                prob += (GB_binary_X[(x,i,h + 2)] - GB_binary_X[(x,i,h + 1)]) - (GB_binary_X[(x,i,h + 1)] - GB_binary_X[(x,i,h)]) >= -1
                prob += (GB_binary_X[(x,i,h + 2)] - GB_binary_X[(x,i,h + 1)]) - (GB_binary_X[(x,i,h + 1)] - GB_binary_X[(x,i,h)]) <= 1
                # 电锅炉
                prob += EB_binary_X[(x,i,h)] * EB_capacity_X[x][i] * EB_minplr_X[x][i] <= EB_hea_output_X[(x,i,h)]
                prob += EB_hea_output_X[(x,i,h)] <= EB_capacity_X[x][i] * EB_binary_X[(x,i,h)]
                prob += EB_hea_output_X[(x,i,h)] == EB_ele_comsuption_X[(x,i,h)] * eff_EB_X[x][i]
                prob += (EB_binary_X[(x,i,h + 2)] - EB_binary_X[(x,i,h + 1)]) - (EB_binary_X[(x,i,h + 1)] - EB_binary_X[(x,i,h)]) >= -1
                prob += (EB_binary_X[(x,i,h + 2)] - EB_binary_X[(x,i,h + 1)]) - (EB_binary_X[(x,i,h + 1)] - EB_binary_X[(x,i,h)]) <= 1
    # 地源热泵
    for x in range(13):
        for h in range(24):
            prob += GSHP_binary_X[(x,0,h)] + GSHP_binary_X[(x,1,h)] + GSHP_binary_X[(x,2,h)] + GSHP_binary_X[(x,3,h)] + GSHP_binary_X[(x,4,h)] + GSHP_binary_X[(x,5,h)] == sum(GSHP_select_configured_X[x])
        prob += sum(GSHP_select_configured_X[x]) <= 6
    # 燃气内燃机
    for x in range(13):
        for h in range(24):
            prob += GICE_binary_X[(x,0,h)] + GICE_binary_X[(x,1,h)] + GICE_binary_X[(x,2,h)] + GICE_binary_X[(x,3,h)] + GICE_binary_X[(x,4,h)] + GICE_binary_X[(x,5,h)] == sum(GICE_select_configured_X[x])
            prob += GICE_ele_output_X[(x,0,h)] + GICE_ele_output_X[(x,1,h)] + GICE_ele_output_X[(x,2,h)] + GICE_ele_output_X[(x,3,h)] + GICE_ele_output_X[(x,4,h)] + GICE_ele_output_X[(x,5,h)] >= GICEEle_sell_grid_X[(x,h)]
        prob += sum(GICE_select_configured_X[x]) <= 6
    # 燃气锅炉
    for x in range(13):
        for h in range(24):
            prob += GB_binary_X[(x,0,h)] + GB_binary_X[(x,1,h)] + GB_binary_X[(x,2,h)] + GB_binary_X[(x,3,h)] + GB_binary_X[(x,4,h)] + GB_binary_X[(x,5,h)] == sum(GB_select_configured_X[x])
        prob += sum(GB_select_configured_X[x]) <= 6
    # 电锅炉
    for x in range(13):
        for h in range(24):
            prob += EB_binary_X[(x,0,h)] + EB_binary_X[(x,1,h)] + EB_binary_X[(x,2,h)] + EB_binary_X[(x,3,h)] + EB_binary_X[(x,4,h)] + EB_binary_X[(x,5,h)] == sum(EB_select_configured_X[x])
        prob += sum(EB_select_configured_X[x]) <= 6
    # 光热
    for x in range(13):
        prob += HC_area_X[x] >= 0
        prob += HC_area_X[x] <= preselection_HC_X[x] * M
        for h in range(24):
            prob += HC_output_X[(x,h)] == HC_area_X[x] * I_solar[h] * eff_HC_X[x]
    # 水蓄能
    for x in range(13):
        prob += WTSH_voulme_X[x] >= 0
        prob += WTSH_voulme_X[x] <= preselection_WTSH_X[x] * M
        prob += WTSH_X[(x,0)] == WTSH_in_X[(x,0)] * eff_Hin_X[x]- WTSH_out_X[(x,0)] * (1/eff_Hout_X[x])
        for h in range(24):
            prob += WTSH_X[(x,h)] <= WTSH_voulme_X[x] * density_H
            prob += WTSH_in_X[(x,h)] <= WTSH_voulme_X[x] * ratio_Hin_X[x]
            prob += WTSH_out_X[(x,h)] <= WTSH_voulme_X[x] * ratio_Hout_X[x]
            if h >= 1:
                prob += WTSH_X[(x,h)] == WTSH_X[(x,h - 1)] * (1 - WTSH_loss_X[x]) + WTSH_in_X[(x,h)] * eff_Hin_X[x] - WTSH_out_X[(x,h)] *(1/ eff_Hout_X[x])
            prob += (WTSH_in_X[(x,h + 2)] - WTSH_in_X[(x,h + 1)]) - (WTSH_in_X[(x,h + 1)] - WTSH_in_X[(x,h)]) >= -1.5 * WTSH_in_X[(x,h + 1)]
            prob += (WTSH_out_X[(x,h + 2)] - WTSH_out_X[(x,h + 1)]) - (WTSH_out_X[(x,h + 1)] - WTSH_out_X[(x,h)]) >= -1.5 * WTSH_out_X[(x,h + 1)]
        prob += sum(WTSH_in_X[(x,h)] for h in range(24)) * eff_Hin_X[x] - sum(WTSH_X[(x,h)] for h in range(24)) * WTSH_loss_X[x] == sum(WTSH_out_X[(x,h)] for h in range(24)) * (1/eff_Hout_X[x])
    # 换热器
    for x in range(13):
        for h in range(24):
            prob += HE_hea_output_X[(x,h)] == input_block[(x,h)] * eff_pipenet_X[x] * eff_HE_X[x]

    # 用户能源站
    for x in range(13):
        for y in range(3):
            # 用户能源站预选约束
            prob += preselection_GSHP_Y[x][y] >= 0
            prob += preselection_GSHP_Y[x][y] <= preselection_user[x][y]
            prob += preselection_GICE_Y[x][y] >= 0
            prob += preselection_GICE_Y[x][y] <= preselection_user[x][y]
            prob += preselection_GB_Y[x][y] >= 0
            prob += preselection_GB_Y[x][y] <= preselection_user[x][y]
            prob += preselection_EB_Y[x][y] >= 0
            prob += preselection_EB_Y[x][y] <= preselection_user[x][y]
            prob += preselection_HC_Y[x][y] >= 0
            prob += preselection_HC_Y[x][y] <= preselection_user[x][y]
            prob += preselection_WTSH_Y[x][y] >= 0
            prob += preselection_WTSH_Y[x][y] <= preselection_user[x][y]
            for i in range(6):
                # 地源热泵
                prob += GSHP_select_configured_Y[x][y][i] >= 0
                prob += GSHP_select_configured_Y[x][y][i] <= preselection_GSHP_Y[x][y]
                # 燃气内燃机
                prob += GICE_select_configured_Y[x][y][i] >= 0
                prob += GICE_select_configured_Y[x][y][i] <= preselection_GICE_Y[x][y]
                # 燃气锅炉
                prob += GB_select_configured_Y[x][y][i] >= 0
                prob += GB_select_configured_Y[x][y][i] <= preselection_GB_Y[x][y]
                # 电锅炉
                prob += EB_select_configured_Y[x][y][i] >= 0
                prob += EB_select_configured_Y[x][y][i] <= preselection_EB_Y[x][y]
                for h in range(24):
                    # 地源热泵
                    prob += GSHP_binary_Y[(x,y,i,h)] * GSHP_capacity_Y[x][y][i] * GSHP_minplr_Y[x][y][i] <= GSHP_output_Y[(x,y,i,h)]
                    prob += GSHP_output_Y[(x,y,i,h)] <= GSHP_capacity_Y[x][y][i] * GSHP_binary_Y[(x,y,i,h)]
                    prob += GSHP_ele_comumption_Y[(x,y,i,h)] == GSHP_output_Y[(x,y,i,h)] * (1/eff_GSHP_Y[x][y][i]) * (1 + GSHP_pump_ele_Y[x][y][i])
                    prob += (GSHP_binary_Y[(x,y,i,h + 2)] - GSHP_binary_Y[(x,y,i,h + 1)]) - (GSHP_binary_Y[(x,y,i,h + 1)] - GSHP_binary_Y[(x,y,i,h)]) >= -1
                    prob += (GSHP_binary_Y[(x,y,i,h + 2)] - GSHP_binary_Y[(x,y,i,h + 1)]) - (GSHP_binary_Y[(x,y,i,h + 1)] - GSHP_binary_Y[(x,y,i,h)]) <= 1
                    # 燃气内燃机
                    prob += GICE_binary_Y[(x,y,i,h)] * GICE_rated_power_Y[x][y][i] * GICE_minplr_Y[x][y][i] <= GICE_ele_output_Y[(x,y,i,h)]
                    prob += GICE_ele_output_Y[(x,y,i,h)] <= GICE_rated_power_Y[x][y][i] * GICE_binary_Y[(x,y,i,h)]
                    prob += GICE_ele_output_Y[(x,y,i,h)] == GICE_gas_comsuption_Y[(x,y,i,h)] * effE_GICE_Y[x][y][i] * q_gas / 3600
                    prob += GICE_hea_output_Y[(x,y,i,h)] == GICE_gas_comsuption_Y[(x,y,i,h)] * effH_GICE_Y[x][y][i] * q_gas / 3600
                    prob += (GICE_binary_Y[(x,y,i,h + 2)] - GICE_binary_Y[(x,y,i,h + 1)]) - (GICE_binary_Y[(x,y,i,h + 1)] - GICE_binary_Y[(x,y,i,h)]) >= -1
                    prob += (GICE_binary_Y[(x,y,i,h + 2)] - GICE_binary_Y[(x,y,i,h + 1)]) - (GICE_binary_Y[(x,y,i,h + 1)] - GICE_binary_Y[(x,y,i,h)]) <= 1
                    # 燃气锅炉
                    prob += GB_binary_Y[(x,y,i,h)] * GB_capacity_Y[x][y][i] * GB_minplr_Y[x][y][i] <= GB_hea_output_Y[(x,y,i,h)]
                    prob += GB_hea_output_Y[(x,y,i,h)] <= GB_capacity_Y[x][y][i] * GB_binary_Y[(x,y,i,h)]
                    prob += GB_hea_output_Y[(x,y,i,h)] == GB_gas_comsuption_Y[(x,y,i,h)] * q_gas * eff_GB_Y[x][y][i] / 3600
                    prob += (GB_binary_Y[(x,y,i,h + 2)] - GB_binary_Y[(x,y,i,h + 1)]) - (GB_binary_Y[(x,y,i,h + 1)] - GB_binary_Y[(x,y,i,h)]) >= -1
                    prob += (GB_binary_Y[(x,y,i,h + 2)] - GB_binary_Y[(x,y,i,h + 1)]) - (GB_binary_Y[(x,y,i,h + 1)] - GB_binary_Y[(x,y,i,h)]) <= 1
                    # 电锅炉
                    prob += EB_binary_Y[(x,y,i,h)] * EB_capacity_Y[x][y][i] * EB_minplr_Y[x][y][i] <= EB_hea_output_Y[(x,y,i,h)]
                    prob += EB_hea_output_Y[(x,y,i,h)] <= EB_capacity_Y[x][y][i] * EB_binary_Y[(x,y,i,h)]
                    prob += EB_hea_output_Y[(x,y,i,h)] == EB_ele_comsuption_Y[(x,y,i,h)] * eff_EB_Y[x][y][i]
                    prob += (EB_binary_Y[(x,y,i,h + 2)] - EB_binary_Y[(x,y,i,h + 1)]) - (EB_binary_Y[(x,y,i,h + 1)] - EB_binary_Y[(x,y,i,h)]) >= -1
                    prob += (EB_binary_Y[(x,y,i,h + 2)] - EB_binary_Y[(x,y,i,h + 1)]) - (EB_binary_Y[(x,y,i,h + 1)] - EB_binary_Y[(x,y,i,h)]) <= 1
    # 地源热泵
    for x in range(13):
        for y in range(3):
            for h in range(24):
                prob += GSHP_binary_Y[(x,y,0,h)] + GSHP_binary_Y[(x,y,1,h)] + GSHP_binary_Y[(x,y,2,h)] + GSHP_binary_Y[(x,y,3,h)] + GSHP_binary_Y[(x,y,4,h)] + GSHP_binary_Y[(x,y,5,h)] == sum(GSHP_select_configured_Y[x][y])
        prob += sum(GSHP_select_configured_Y[x][y]) <= 6
    # 燃气内燃机
    for x in range(13):
        for y in range(3):
            for h in range(24):
                prob += GICE_binary_Y[(x,y,0,h)]  + GICE_binary_Y[(x,y,1,h)] + GICE_binary_Y[(x,y,2,h)] + GICE_binary_Y[(x,y,3,h)] +GICE_binary_Y[(x,y,4,h)] + GICE_binary_Y[(x,y,5,h)] == sum(GICE_select_configured_Y[x][y])
                prob += GICE_ele_output_Y[(x,y,0,h)] + GICE_ele_output_Y[(x,y,1,h)] + GICE_ele_output_Y[(x,y,2,h)] + GICE_ele_output_Y[(x,y,3,h)] + GICE_ele_output_Y[(x,y,4,h)] + GICE_ele_output_Y[(x,y,5,h)] >= GICEEle_sell_grid_Y[(x,y,h)]
        prob += sum(GICE_select_configured_Y[x][y]) <= 6
    # 燃气锅炉
    for x in range(13):
        for y in range(3):
            for h in range(24):
                prob += GB_binary_Y[(x,y,0,h)]  + GB_binary_Y[(x,y,1,h)] + GB_binary_Y[(x,y,2,h)] + GB_binary_Y[(x,y,3,h)] + GB_binary_Y[(x,y,4,h)] + GB_binary_Y[(x,y,5,h)] == sum(GB_select_configured_Y[x][y])
        prob += sum(GB_select_configured_Y[x][y]) <= 6
    # 电锅炉
    for x in range(13):
        for y in range(3):
            for h in range(24):
                prob += EB_binary_Y[(x,y,0,h)] + EB_binary_Y[(x,y,1,h)] + EB_binary_Y[(x,y,2,h)] + EB_binary_Y[(x,y,3,h)] + EB_binary_Y[(x,y,4,h)] + EB_binary_Y[(x,y,5,h)] == sum(EB_select_configured_Y[x][y])
        prob += sum(EB_select_configured_Y[x][y]) <= 6
    # 光热
    for x in range(13):
        for y in range(3):
            prob += HC_area_Y[x][y] >= 0
            prob += HC_area_Y[x][y] <= preselection_HC_Y[x][y] * M
            for h in range(24):
                prob += HC_output_Y[(x,y,h)] == HC_area_Y[x][y] * I_solar[h] * eff_HC_Y[x][y]
    # 水蓄能
    for x in range(13):
        for y in range(3):
            prob += WTSH_voulme_Y[x][y] >= 0
            prob += WTSH_voulme_Y[x][y] <= preselection_WTSH_Y[x][y] * M
            prob += WTSH_Y[(x,y,0)] == WTSH_in_Y[(x,y,0)] * eff_Hin_Y[x][y] - WTSH_out_Y[(x,y,0)] * (1/eff_Hout_Y[x][y])
            for h in range(24):
                prob += WTSH_Y[(x,y,h)] <= WTSH_voulme_Y[x][y] * density_H
                prob += WTSH_in_Y[(x,y,h)] <= WTSH_voulme_Y[x][y] * ratio_Hin_Y[x][y]
                prob += WTSH_out_Y[(x,y,h)] <= WTSH_voulme_Y[x][y] * ratio_Hout_Y[x][y]
                if h >= 1:
                    prob += WTSH_Y[(x,y,h)] == WTSH_Y[(x,y,h - 1)] * (1 - WTSH_loss_Y[x][y]) + WTSH_in_Y[(x,y,h)] * eff_Hin_Y[x][y] - WTSH_out_Y[(x,y,h)] * (1/eff_Hout_Y[x][y])
                prob += (WTSH_in_Y[(x,y,h + 2)] - WTSH_in_Y[(x,y,h + 1)]) - (WTSH_in_Y[(x,y,h + 1)]- WTSH_in_Y[(x,y,h)]) >= -1.5 * WTSH_in_Y[(x,y,h + 1)]
                prob += (WTSH_out_Y[(x,y,h + 2)] - WTSH_out_Y[(x,y,h + 1)]) - (WTSH_out_Y[(x,y,h + 1)] - WTSH_out_Y[(x,y,h)]) >= -1.5 * WTSH_out_Y[(x,y,h + 1)]
            prob += sum(WTSH_in_Y[(x,y,h)] for h in range(24)) * eff_Hin_Y[x][y] - sum(WTSH_Y[(x,y,h)] for h in range(24)) * WTSH_loss_Y[x][y] == sum(WTSH_out_Y[(x,y,h)] for h in range(24)) * (1/eff_Hout_Y[x][y])

    # 结果展示用
    # 区域能源站总购电费用、总售电费用、总购燃气费用、日总耗电量、日总售电量、日总耗燃气量、小时总耗电量、小时总售电量, 小时总耗燃气量
    Ele_cost_regi = LpVariable('Ele_cost_regi')
    Ele_sell_regi = LpVariable('Ele_sell_regi')
    Gas_cost_regi = LpVariable('Gas_cost_regi')
    Ele_daily_regi = LpVariable('Ele_daily_regi')
    Ele_sell_daily_regi = LpVariable('Ele_sell_daily_regi')
    Gas_daily_regi = LpVariable('Gas_daily_regi')
    Ele_hour_regi = LpVariable.dicts('Ele_hour_regi', (h for h in hh), lowBound=0)
    Ele_sell_hour_regi = LpVariable.dicts('Ele_sell_hour_regi', (h for h in hh), lowBound=0)
    Gas_hour_regi = LpVariable.dicts('Gas_hour_regi', (h for h in hh), lowBound=0)
    # 街区能源站总购电费用、总售电费用、总购燃气费用、日总耗电量、日总售电量、日总耗燃气量、小时总耗电量、小时总售电量, 小时总耗燃气量
    # Ele_hour_block(x, h), Ele_sell_hour_block(x, h), Gas_hour_block(x, h);
    Ele_cost_block = LpVariable.dicts('Ele_cost_block', (x for x in xx), lowBound=0)
    Ele_sell_block = LpVariable.dicts('Ele_sell_block', (x for x in xx), lowBound=0)
    Gas_cost_block = LpVariable.dicts('Gas_cost_block', (x for x in xx), lowBound=0)
    Ele_daily_block = LpVariable.dicts('Ele_daily_block', (x for x in xx), lowBound=0)
    Ele_sell_daily_block = LpVariable.dicts('Ele_sell_daily_block', (x for x in xx), lowBound=0)
    Gas_daily_block = LpVariable.dicts('Gas_daily_block', (x for x in xx), lowBound=0)
    Ele_hour_block = LpVariable.dicts('Ele_hour_block', ((x, h) for x in xx for h in hh), lowBound=0)
    Ele_sell_hour_block = LpVariable.dicts('Ele_sell_hour_block', ((x, h) for x in xx for h in hh), lowBound=0)
    Gas_hour_block = LpVariable.dicts('Gas_hour_block', ((x, h) for x in xx for h in hh), lowBound=0)
    # 用户能源站总购电费用、总售电费用、总购燃气费用、日总耗电量、日总售电量、日总耗燃气量、小时总耗电量、小时总售电量, 小时总耗燃气量
    Ele_cost_user = LpVariable.dicts('Ele_cost_user', ((x, y) for x in xx for y in yy), lowBound=0)
    Ele_sell_user = LpVariable.dicts('Ele_sell_user', ((x, y) for x in xx for y in yy), lowBound=0)
    Gas_cost_user = LpVariable.dicts('Gas_cost_user', ((x, y) for x in xx for y in yy), lowBound=0)
    Ele_daily_user = LpVariable.dicts('Ele_daily_user', ((x, y) for x in xx for y in yy), lowBound=0)
    Ele_sell_daily_user = LpVariable.dicts('Ele_sell_daily_user', ((x, y) for x in xx for y in yy), lowBound=0)
    Gas_daily_user = LpVariable.dicts('Gas_daily_user', ((x, y) for x in xx for y in yy), lowBound=0)
    Ele_hour_user = LpVariable.dicts('Ele_hour_user', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)
    Ele_sell_hour_user = LpVariable.dicts('Ele_sell_hour_user', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)
    Gas_hour_user = LpVariable.dicts('Gas_hour_user', ((x, y, h) for x in xx for y in yy for h in hh), lowBound=0)

    # 区域能源站
    prob += Ele_cost_regi == sum(Ele_buy_grid[h] * Ele_price[h] for h in hh)
    prob += Ele_sell_regi == sum(GICEEle_sell_grid[h] * Ele_sell_price for h in hh)
    prob += Gas_cost_regi == sum(GICE_gas_comsuption[(i,h)] + GB_gas_comsuption[(i,h)] * Gas_price for i in ii for h in hh)
    prob += Ele_daily_regi == sum(Ele_buy_grid[h] for h in hh)
    prob += Ele_sell_daily_regi == sum(GICEEle_sell_grid[h] for h in hh)
    prob += Gas_daily_regi == sum(GICE_gas_comsuption[(i,h)] + GB_gas_comsuption[(i,h)] for i in ii for h in hh)
    for h in hh:
        prob += Ele_hour_regi[h] == Ele_buy_grid[h]
        prob += Ele_sell_hour_regi[h] == GICETotalEle_sell_grid[h]
        prob += Gas_hour_regi[h] == sum(GICE_gas_comsuption[(i, h)] + GB_gas_comsuption[(i, h)] for i in ii)
        # 区域能源站输给街区能源站热量 / 区域能源站的小时热负荷
        prob += output_regi[h] == sum(input_block[(x,h)] for x in xx)
        # 区域能源站输给街区能源站电量
        prob += output_ele_regi[h] == sum(input_ele_block[(x,h)] for x in xx)

    # 街区能源站
    for x in xx:
        prob += Ele_cost_block[x] == sum(Ele_buy_grid_X[(x, h)] * Ele_price[h] for h in hh)
        prob += Ele_sell_block[x] == sum(GICEEle_sell_grid_X[(x, h)] * Ele_sell_price for h in hh)
        prob += Gas_cost_block[x] == sum(GICE_gas_comsuption_X[(x, i, h)] + GB_gas_comsuption_X[(x, i, h)] * Gas_price for i in ii for h in hh)
        prob += Ele_daily_block[x] == sum(Ele_buy_grid_X[(x,h)] for h in hh)
        prob += Ele_sell_daily_block[x] == sum(GICEEle_sell_grid_X[(x,h)] for h in hh)
        prob += Gas_daily_block[x] == sum(GICE_gas_comsuption_X[(x,i,h)] + GB_gas_comsuption_X[(x,i,h)] for i in ii for h in hh)
    for x in xx:
        for h in hh:
            prob += Ele_hour_block[(x,h)] == Ele_buy_grid_X[(x,h)]
            prob += Ele_sell_hour_block[(x,h)] == GICEEle_sell_grid_X[(x,h)]
            prob += Gas_hour_block[(x,h)] == sum(GICE_gas_comsuption_X[(x,i,h)] + GB_gas_comsuption_X[(x, i, h)] for i in ii)
            # 每个街区能源站输给用户能源站总热量 / 每个街区能源站的小时负荷
            prob += output_block[(x,h)] == sum(input_user[(x,y,h)] for y in yy)
            # 街区能源站输给用户能源站底电量
            prob += output_ele_block[(x, h)] == sum(input_ele_user[(x, y, h)] for y in yy)

    # 用户能源站
    for x in xx:
        for y in yy:
            prob += Ele_cost_user[(x, y)] == sum(Ele_buy_grid_Y[(x, y, h)] * Ele_price[h] for h in hh)
            prob += Ele_sell_user[(x, y)] == sum(GICEEle_sell_grid_Y[(x, y, h)] * Ele_sell_price for h in hh)
            prob += Gas_cost_user[(x, y)] == sum(GICE_gas_comsuption_Y[(x, y,i, h)] + GB_gas_comsuption_Y[(x, y, i,h)] * Gas_price for i in ii for h in hh)
            prob += Ele_daily_user[(x, y)] == sum(Ele_buy_grid_Y[(x, y, h)] for h in hh)
            prob += Ele_sell_daily_user[(x, y)] == sum(GICEEle_sell_grid_Y[(x, y, h)] for h in hh)
            prob += Gas_daily_user[(x, y)] == sum(GICE_gas_comsuption_Y[(x, y,i, h)] + GB_gas_comsuption_Y[(x, y, i, h)] for i in ii for h in hh)
    for x in xx:
        for y in yy:
            for h in hh:
                prob += Ele_hour_user[(x, y, h)] == Ele_buy_grid_Y[(x, y, h)]
                prob += Ele_sell_hour_user[(x, y, h)] == GICEEle_sell_grid_Y[(x, y, h)]
                prob += Gas_hour_user[(x, y, h)] == sum(GICE_gas_comsuption_Y[(x, y, i, h)] + GB_gas_comsuption_Y[(x, y, i,h)] for i in ii)
                # 用户能源站从街区能源站获得的热量
                prob += input_user_Y[(x, y, h)] == input_user[(x, y, h)]* eff_pipenet_Y[(x, y)]

    # 求解
    prob.solve(solver=CPLEX())
    # 查看解的状态
    print("Status:", LpStatus[prob.status])
    # 查看解
    # for v in prob.variables():
    #     print(v.name, "=", v.varValue)
    print(value(prob.objective))


def yh(clid, date, cityid, zxfs):
    # 读取配置文件
    cp = configparser.ConfigParser()
    cp.read("config.ini")
    user = cp.get("db", "user")
    password = cp.get("db", "password")
    database = cp.get("db", "database")
    host = cp.get("db", "host")
    port = cp.get("db", "port")
    # 对应的数据库
    conn = mysql.connector.connect(host=host, port=port, user=user, password=password, database=database)
    cursor = conn.cursor()

    date = datetime.datetime.strptime(date, "%Y-%m-%d")
    now = datetime.datetime.strptime(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '%Y-%m-%d %H:%M:%S')
    # 根据策略id查询类别和是否启动
    # 生效时间 <= 执行当前时间, 取最大生效时间那条数据的version
    cursor.execute(
        'select lsxgid,clid,version,lb,IFNULL(sfqy,1),yhmb,sdjg,rqjg,rqrz from taskinvoke_yhdd_clzdls where effecttime = (select max(effecttime) from taskinvoke_yhdd_clzdls where clid = %s and effecttime <= %s)',
        (clid, now))
    nyz = cursor.fetchall()[0]
    # print(nyz)
    # 判断是集中式系统还是分布式系统 1集中式 2分布式
    if nyz[3] == '1':
        # 集中式系统
        # 判断是否启用 1启用 0未启用
        if nyz[4] == '1':
            # 启用
            preselection_regi = 1  # 区域能源站预选
            # 判断优化目标  01位费用  02为排放
            if nyz[5] == '01':
                # 费用为优化目标

                # 找最顶层能源站（区域能源站）
                cursor.execute('select nyzid from taskinvoke_yhdd_cynyz where clid = %s and pid is null', (clid,))
                qyid = cursor.fetchall()[0][0]  # 区域能源站
                # 找中层能源站（街区能源站）
                cursor.execute('select nyzid from taskinvoke_yhdd_cynyz where clid = %s and pid = %s', (clid, qyid))
                jqid = cursor.fetchall()
                # 找最底层能源站（用户能源站）
                jqids = []  # 街区能源站
                yhids = []  # 用户能源站
                for jq in jqid:
                    jqids.append(jq[0])
                    cursor.execute('select nyzid from taskinvoke_yhdd_cynyz where clid = %s and pid = %s',
                                   (clid, jq[0]))
                    yhid = cursor.fetchall()
                    yhs = []
                    for yh in yhid:
                        yhs.append(yh[0])
                    yhids.append(yhs)
                if len(jqids) != 13:
                    num = 13 - len(jqids)
                    for n in range(num):
                        jqids.append('0')
                if len(yhids) != 13:
                    for u in yhids:
                        if len(u) != 3:
                            cha = 3 - len(u)
                            u.append('0')
                    a = 13 - len(yhids)
                    b = ['0', '0', '0']
                    for c in range(a):
                        yhids.append(b)
                # print(qyid)
                # print(jqids)
                # print(yhids)
                # print(len(yhids))

                # 查询热负荷
                Load_hea = []  # 热负荷
                for yid in yhids:
                    y = []
                    for yi in yid:
                        fhycid = ''.join(str(date)[0:10].split('-')) + yi
                        cursor.execute('select fh from taskinvoke_fhyc_dq where fhycid = %s', (fhycid,))
                        fhs = cursor.fetchall()
                        # print(fhs)
                        list = []
                        if len(fhs) == 0:
                            # a0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                            for w in range(24):
                                list.append(0)
                        for f in fhs:
                            list.append(float(f[0]))
                        y.append(list)
                    Load_hea.append(y)
                # 查询辐射强度
                Ndate = date.timetuple().tm_yday
                I_solar = []  # 太阳辐射
                for i in range(0, 24):
                    Ih = getIhdata(Ndate, i, date, cityid, cursor)
                    I_solar.append(Ih)
                # 查询24小时购电价格
                cursor.execute('select sjd,dj from t_bm_fsdj where sfid = %s order by sjd', (cityid[0:2],))
                gdj = cursor.fetchall()
                Ele_price = []  # 24小时购电价格
                for gd in gdj:
                    Ele_price.append(float(gd[1]))
                # 查询售电价格
                Ele_sell_price = float(nyz[6])
                # 查询燃气价格
                Gas_price = float(nyz[7])
                # 查询燃气热值
                q_gas = nyz[8]

                # 查询参与设备
                # 区域能源站
                # 热源预选
                preselection_GSHP = getYx(nyz[0], qyid, '01', cursor)  # 地源热泵预选
                preselection_GICE = getYx(nyz[0], qyid, '02', cursor)  # 燃气内燃机预选
                preselection_GB = getYx(nyz[0], qyid, '03', cursor)  # 燃气锅炉预选
                preselection_EB = getYx(nyz[0], qyid, '04', cursor)  # 电锅炉预选
                # 地源热泵
                GSHP_select_configured = getTs(nyz[0], qyid, '01', cursor)  # 地源热泵台数
                qydyrbcs = getCs('sbrl', 'sbxl', 'zdfhl', 'sbhdxs', nyz[0], qyid, '01', cursor, '1')
                GSHP_capacity = qydyrbcs[0]  # 地源热泵容量
                eff_GSHP = qydyrbcs[1]  # 地源热泵效率
                GSHP_minplr = qydyrbcs[2]  # 地源热泵最低负荷率
                GSHP_pump_ele = qydyrbcs[3]  # 地源热泵水泵耗电系数
                # 燃气内燃机
                GICE_select_configured = getTs(nyz[0], qyid, '02', cursor)  # 燃气内燃机台数
                qyrqnrjcs = getCs('edfdl', 'fdxl', 'zrxl', 'zdfhl', nyz[0], qyid, '02', cursor, '1')
                GICE_rated_power = qyrqnrjcs[0]  # 燃气内燃机额定发电量
                effE_GICE = qyrqnrjcs[1]  # 燃气内燃机发电效率
                effH_GICE = qyrqnrjcs[2]  # 燃气内燃机制热效率
                GICE_minplr = qyrqnrjcs[3]  # 燃气内燃机最低负荷率
                # 燃气锅炉
                GB_select_configured = getTs(nyz[0], qyid, '03', cursor)  # 燃气锅炉台数
                qyrqglcs = getCs2('sbrl', 'sbxl', 'zdfhl', nyz[0], qyid, '03', cursor, '1')
                GB_capacity = qyrqglcs[0]  # 燃气锅炉容量
                eff_GB = qyrqglcs[1]  # 燃气锅炉效率
                GB_minplr = qyrqglcs[2]  # 燃气锅炉最低负荷率
                # 电锅炉
                EB_select_configured = getTs(nyz[0], qyid, '04', cursor)  # 电锅炉台数
                qydglcs = getCs2('sbrl', 'sbxl', 'zdfhl', nyz[0], qyid, '04', cursor, '1')
                EB_capacity = qydglcs[0]  # 电锅炉容量
                eff_EB = qydglcs[1]  # 电锅炉效率
                EB_minplr = qydglcs[2]  # 电锅炉最低负荷率

                # 街区能源站
                preselection_block = []  # 街区能源站预选
                for jq in jqids:
                    if jq != '0':
                        preselection_block.append(1)
                    else:
                        preselection_block.append(0)
                # 热源预选
                preselection_GSHP_X = []  # 地源热泵预选
                preselection_GICE_X = []  # 燃气内燃机预选
                preselection_GB_X = []  # 燃气锅炉预选
                preselection_EB_X = []  # 电锅炉预选
                preselection_HC_X = []  # 光热预选
                preselection_WTSH_X = []  # 水蓄能预选 preselection_block, preselection_GSHP_X, preselection_GICE_X, preselection_GB_X, preselection_EB_X, preselection_HC_X, preselection_WTSH_X
                eff_pipenet_X = []  # 管网输配效率
                cursor.execute(
                    'select gwspxl from t_yhdd_gwspxl where id = 1')
                spxl = cursor.fetchall()[0][0]
                for jqs in jqids:
                    preselection_GSHP_X.append(getYx(nyz[0], jqs, '01', cursor))
                    preselection_GICE_X.append(getYx(nyz[0], jqs, '02', cursor))
                    preselection_GB_X.append(getYx(nyz[0], jqs, '03', cursor))
                    preselection_EB_X.append(getYx(nyz[0], jqs, '04', cursor))
                    preselection_HC_X.append(getYx(nyz[0], jqs, '05', cursor))
                    preselection_WTSH_X.append(getYx(nyz[0], jqs, '06', cursor))
                    if jqs != '0':
                        eff_pipenet_X.append(float(spxl))
                    else:
                        eff_pipenet_X.append(0)
                GSHP_select_configured_X = []  # 地源热泵台数
                GSHP_capacity_X = []  # 地源热泵容量
                eff_GSHP_X = []  # 地源热泵效率
                GSHP_minplr_X = []  # 地源热泵最低负荷率
                GSHP_pump_ele_X = []  # 地源热泵水泵耗电系数 GSHP_select_configured_X, GSHP_capacity_X, eff_GSHP_X, GSHP_minplr_X, GSHP_pump_ele_X
                GICE_select_configured_X = []  # 燃气内燃机台数
                GICE_rated_power_X = []  # 燃气内燃机额定发电量
                effE_GICE_X = []  # 燃气内燃机发电效率
                effH_GICE_X = []  # 燃气内燃机制热效率
                GICE_minplr_X = []  # 燃气内燃机最低负荷率 GICE_select_configured_X, GICE_rated_power_X, effE_GICE_X, effH_GICE_X, GICE_minplr_X
                GB_select_configured_X = []  # 燃气锅炉台数
                GB_capacity_X = []  # 燃气锅炉容量
                eff_GB_X = []  # 燃气锅炉效率
                GB_minplr_X = []  # 燃气锅炉最低负荷率 GB_select_configured_X, GB_capacity_X, eff_GB_X, GB_minplr_X
                EB_select_configured_X = []  # 电锅炉台数
                EB_capacity_X = []  # 电锅炉容量
                eff_EB_X = []  # 电锅炉效率
                EB_minplr_X = []  # 电锅炉最低负荷率 EB_select_configured_X, EB_capacity_X, eff_EB_X, EB_minplr_X
                HC_area_X = []  # 光热面积
                eff_HC_X = []  # 光热效率 HC_area_X, eff_HC_X
                WTSH_voulme_X = []  # 水蓄热容量
                eff_Hin_X = []  # 水蓄热储能效率
                eff_Hout_X = []  # 水蓄热释能效率
                WTSH_loss_X = []  # 水蓄热自损率
                ratio_Hin_X = []  # 水蓄热储能速率
                ratio_Hout_X = []  # 水蓄热释能速率 WTSH_voulme_X, eff_Hin_X, eff_Hout_X, WTSH_loss_X, ratio_Hin_X, ratio_Hout_X
                HE_capacity_X = []  # 换热器容量
                eff_HE_X = []  # 换热器效率 eff_pipenet_X, eff_HE_X
                for jqs in jqids:
                    # 地源热泵
                    ts1 = getTs(nyz[0], jqs, '01', cursor)
                    GSHP_select_configured_X.append(ts1)
                    cs1 = getCs('sbrl', 'sbxl', 'zdfhl', 'sbhdxs', nyz[0], jqs, '01', cursor, '2')
                    rl1 = cs1[0]
                    GSHP_capacity_X.append(rl1)
                    xl1 = cs1[1]
                    eff_GSHP_X.append(xl1)
                    zdfhl1 = cs1[2]
                    GSHP_minplr_X.append(zdfhl1)
                    rbsbhdxs1 = cs1[3]
                    GSHP_pump_ele_X.append(rbsbhdxs1)
                    # 燃气内燃机
                    ts2 = getTs(nyz[0], jqs, '02', cursor)
                    GICE_select_configured_X.append(ts2)
                    cs2 = getCs('edfdl', 'fdxl', 'zrxl', 'zdfhl', nyz[0], jqs, '02', cursor, '2')
                    dl2 = cs2[0]
                    GICE_rated_power_X.append(dl2)
                    fdxl2 = cs2[1]
                    effE_GICE_X.append(fdxl2)
                    zrxl2 = cs2[2]
                    effH_GICE_X.append(zrxl2)
                    zdfhl2 = cs2[3]
                    GICE_minplr_X.append(zdfhl2)
                    # 燃气锅炉
                    ts3 = getTs(nyz[0], jqs, '03', cursor)
                    GB_select_configured_X.append(ts3)
                    cs3 = getCs2('sbrl', 'sbxl', 'zdfhl', nyz[0], jqs, '03', cursor, '2')
                    rl3 = cs3[0]
                    GB_capacity_X.append(rl3)
                    xl3 = cs3[1]
                    eff_GB_X.append(xl3)
                    zdfhl3 = cs3[2]
                    GB_minplr_X.append(zdfhl3)
                    # 电锅炉
                    ts4 = getTs(nyz[0], jqs, '04', cursor)
                    EB_select_configured_X.append(ts4)
                    cs4 = getCs2('sbrl', 'sbxl', 'zdfhl', nyz[0], jqs, '04', cursor, '2')
                    rl4 = cs4[0]
                    EB_capacity_X.append(rl4)
                    xl4 = cs4[1]
                    eff_EB_X.append(xl4)
                    zdfhl4 = cs4[2]
                    EB_minplr_X.append(zdfhl4)
                    # 光热
                    cursor.execute(
                        'select grmj,sbxl from taskinvoke_yhdd_clcysbls where lsxgid = %s and nyzid = %s and sblxid = %s',
                        (nyz[0], jqs, '05'))
                    gr = cursor.fetchall()
                    if len(gr) > 0:
                        HC_area_X.append(float(gr[0][0]))
                        eff_HC_X.append(float(gr[0][1]))
                    else:
                        HC_area_X.append(0)
                        eff_HC_X.append(0.45)
                    # 水蓄热
                    cursor.execute(
                        'select sbrl,cnxl,snxl,zsl,cnsl,snsl from taskinvoke_yhdd_clcysbls where lsxgid = %s and nyzid = %s and sblxid = %s',
                        (nyz[0], jqs, '06'))
                    sxr = cursor.fetchall()
                    if len(sxr) > 0:
                        WTSH_voulme_X.append(float(sxr[0][0]))
                        eff_Hin_X.append(float(sxr[0][1]))
                        eff_Hout_X.append(float(sxr[0][2]))
                        WTSH_loss_X.append(float(sxr[0][3]))
                        ratio_Hin_X.append(float(sxr[0][4]))
                        ratio_Hout_X.append(float(sxr[0][5]))
                    else:
                        WTSH_voulme_X.append(0)
                        eff_Hin_X.append(0.87)
                        eff_Hout_X.append(0.87)
                        WTSH_loss_X.append(0.02)
                        ratio_Hin_X.append(0.07)
                        ratio_Hout_X.append(0.07)
                    # 换热器
                    cursor.execute(
                        'select sbrl,sbxl from taskinvoke_yhdd_clcysbls where lsxgid = %s and nyzid = %s and sblxid = %s',
                        (nyz[0], jqs, '07'))
                    hrq = cursor.fetchall()
                    if len(hrq) > 0:
                        HE_capacity_X.append(float(hrq[0][0]))
                        eff_HE_X.append(float(hrq[0][1]))
                    else:
                        HE_capacity_X.append(0)
                        eff_HE_X.append(0)

                # 用户能源站
                preselection_user = []  # 用户能源站预选
                for yhidd in yhids:
                    yyyid = []
                    for yid in yhidd:
                        if yid != '0':
                            yyyid.append(1)
                        else:
                            yyyid.append(0)
                    preselection_user.append(yyyid)
                # 热源预选
                preselection_GSHP_Y = []  # 地源热泵预选
                preselection_GICE_Y = []  # 燃气内燃机预选
                preselection_GB_Y = []  # 燃气锅炉预选
                preselection_EB_Y = []  # 电锅炉预选
                preselection_HC_Y = []  # 光热预选
                preselection_WTSH_Y = []  # 水蓄能预选 preselection_user, preselection_GSHP_Y, preselection_GICE_Y, preselection_GB_Y, preselection_EB_Y, preselection_HC_Y, preselection_WTSH_Y
                eff_pipenet_Y = []  # 管网输配效率
                cursor.execute(
                    'select gwspxl from t_yhdd_gwspxl where id = 1')
                spxl = cursor.fetchall()[0][0]
                for yhidd in yhids:
                    yyhdyrbyx = []
                    yyhrqnrjyx = []
                    yyhrqglyx = []
                    yyhdglyx = []
                    yyhgryx = []
                    yyhsxnyx = []
                    yyhgwspxl = []
                    for yid in yhidd:
                        yyhdyrbyx.append(getYx(nyz[0], yid, '01', cursor))
                        yyhrqnrjyx.append(getYx(nyz[0], yid, '02', cursor))
                        yyhrqglyx.append(getYx(nyz[0], yid, '03', cursor))
                        yyhdglyx.append(getYx(nyz[0], yid, '04', cursor))
                        yyhgryx.append(getYx(nyz[0], yid, '05', cursor))
                        yyhsxnyx.append(getYx(nyz[0], yid, '06', cursor))
                        yyhgwspxl.append(float(spxl))
                    preselection_GSHP_Y.append(yyhdyrbyx)
                    preselection_GICE_Y.append(yyhrqnrjyx)
                    preselection_GB_Y.append(yyhrqglyx)
                    preselection_EB_Y.append(yyhdglyx)
                    preselection_HC_Y.append(yyhgryx)
                    preselection_WTSH_Y.append(yyhsxnyx)
                    eff_pipenet_Y.append(yyhgwspxl)
                GSHP_select_configured_Y = []  # 地源热泵台数
                GSHP_capacity_Y = []  # 地源热泵容量
                eff_GSHP_Y = []  # 地源热泵效率
                GSHP_minplr_Y = []  # 地源热泵最低负荷率
                GSHP_pump_ele_Y = []  # 地源热泵水泵耗电系数 GSHP_select_configured_Y, GSHP_capacity_Y, eff_GSHP_Y, GSHP_minplr_Y, GSHP_pump_ele_Y
                GICE_select_configured_Y = []  # 燃气内燃机台数
                GICE_rated_power_Y = []  # 燃气内燃机额定发电量
                effE_GICE_Y = []  # 燃气内燃机发电效率
                effH_GICE_Y = []  # 燃气内燃机制热效率
                GICE_minplr_Y = []  # 燃气内燃机最低负荷率 GICE_select_configured_Y, GICE_rated_power_Y, effE_GICE_Y, effH_GICE_Y, GICE_minplr_Y
                GB_select_configured_Y = []  # 燃气锅炉台数
                GB_capacity_Y = []  # 燃气锅炉容量
                eff_GB_Y = []  # 燃气锅炉效率
                GB_minplr_Y = []  # 燃气锅炉最低负荷率 GB_select_configured_Y, GB_capacity_Y, eff_GB_Y, GB_minplr_Y
                EB_select_configured_Y = []  # 电锅炉台数
                EB_capacity_Y = []  # 电锅炉容量
                eff_EB_Y = []  # 电锅炉效率
                EB_minplr_Y = []  # 电锅炉最低负荷率 EB_select_configured_Y, EB_capacity_Y, eff_EB_Y, EB_minplr_Y
                HC_area_Y = []  # 光热面积
                eff_HC_Y = []  # 光热效率 HC_area_Y, eff_HC_Y
                WTSH_voulme_Y = []  # 水蓄热容量
                eff_Hin_Y = []  # 水蓄热储能效率
                eff_Hout_Y = []  # 水蓄热释能效率
                WTSH_loss_Y = []  # 水蓄热自损率
                ratio_Hin_Y = []  # 水蓄热储能速率
                ratio_Hout_Y = []  # 水蓄热释能速率 WTSH_voulme_Y, eff_Hin_Y, eff_Hout_Y, WTSH_loss_Y, ratio_Hin_Y, ratio_Hout_Y
                for yhidd in yhids:
                    yyhdyrbts = []
                    yyhdyrbrl = []
                    yyhdyrbxl = []
                    yyhdyrbzdfhl = []
                    yyhdyrbsbhdxs = []
                    yyhrqnrjts = []
                    yyhrqnrjfdl = []
                    yyhrqnrjfdxl = []
                    yyhrqnrjzrxl = []
                    yyhrqnrjzdfhl = []
                    yyhrqglts = []
                    yyhrqglrl = []
                    yyhrqglxl = []
                    yyhrqglzdfhl = []
                    yyhdglts = []
                    yyhdglrl = []
                    yyhdglxl = []
                    yyhdglzdfhl = []
                    yyhgrmj = []
                    yyhgrxl = []
                    yyhsxrrl = []
                    yyhsxrcnxl = []
                    yyhsxrsnxl = []
                    yyhsxrzsl = []
                    yyhsxrcnsl = []
                    yyhsxrsnsl = []
                    for yid in yhidd:
                        # 地源热泵
                        ts1 = getTs(nyz[0], yid, '01', cursor)
                        yyhdyrbts.append(ts1)
                        cs1 = getCs('sbrl', 'sbxl', 'zdfhl', 'sbhdxs', nyz[0], yid, '01', cursor, '3')
                        rl1 = cs1[0]
                        yyhdyrbrl.append(rl1)
                        xl1 = cs1[1]
                        yyhdyrbxl.append(xl1)
                        zdfhl1 = cs1[2]
                        yyhdyrbzdfhl.append(zdfhl1)
                        rbsbhdxs1 = cs1[3]
                        yyhdyrbsbhdxs.append(rbsbhdxs1)
                        # 燃气内燃机
                        ts2 = getTs(nyz[0], yid, '02', cursor)
                        yyhrqnrjts.append(ts2)
                        cs2 = getCs('edfdl', 'fdxl', 'zrxl', 'zdfhl', nyz[0], yid, '02', cursor, '3')
                        dl2 = cs2[0]
                        yyhrqnrjfdl.append(dl2)
                        fdxl2 = cs2[1]
                        yyhrqnrjfdxl.append(fdxl2)
                        zrxl2 = cs2[2]
                        yyhrqnrjzrxl.append(zrxl2)
                        zdfhl2 = cs2[3]
                        yyhrqnrjzdfhl.append(zdfhl2)
                        # 燃气锅炉
                        ts3 = getTs(nyz[0], yid, '03', cursor)
                        yyhrqglts.append(ts3)
                        cs3 = getCs2('sbrl', 'sbxl', 'zdfhl', nyz[0], yid, '03', cursor, '3')
                        rl3 = cs3[0]
                        yyhrqglrl.append(rl3)
                        xl3 = cs3[1]
                        yyhrqglxl.append(xl3)
                        zdfhl3 = cs3[2]
                        yyhrqglzdfhl.append(zdfhl3)
                        # 电锅炉
                        ts4 = getTs(nyz[0], yid, '04', cursor)
                        yyhdglts.append(ts4)
                        cs4 = getCs2('sbrl', 'sbxl', 'zdfhl', nyz[0], yid, '04', cursor, '3')
                        rl4 = cs4[0]
                        yyhdglrl.append(rl4)
                        xl4 = cs4[1]
                        yyhdglxl.append(xl4)
                        zdfhl4 = cs4[2]
                        yyhdglzdfhl.append(zdfhl4)
                        # 光热
                        cursor.execute(
                            'select grmj,sbxl from taskinvoke_yhdd_clcysbls where lsxgid = %s and nyzid = %s and sblxid = %s',
                            (nyz[0], yid, '05'))
                        gr = cursor.fetchall()
                        if len(gr) > 0:
                            yyhgrmj.append(float(gr[0][0]))
                            yyhgrxl.append(float(gr[0][1]))
                        else:
                            yyhgrmj.append(0)
                            yyhgrxl.append(0)
                        # 水蓄热
                        cursor.execute(
                            'select sbrl,cnxl,snxl,zsl,cnsl,snsl from taskinvoke_yhdd_clcysbls where lsxgid = %s and nyzid = %s and sblxid = %s',
                            (nyz[0], yid, '06'))
                        sxr = cursor.fetchall()
                        if len(sxr) > 0:
                            yyhsxrrl.append(float(sxr[0][0]))
                            yyhsxrcnxl.append(float(sxr[0][1]))
                            yyhsxrsnxl.append(float(sxr[0][2]))
                            yyhsxrzsl.append(float(sxr[0][3]))
                            yyhsxrcnsl.append(float(sxr[0][4]))
                            yyhsxrsnsl.append(float(sxr[0][5]))
                        else:
                            yyhsxrrl.append(0)
                            yyhsxrcnxl.append(0.87)
                            yyhsxrsnxl.append(0.87)
                            yyhsxrzsl.append(0.02)
                            yyhsxrcnsl.append(0.09)
                            yyhsxrsnsl.append(0.09)
                    GSHP_select_configured_Y.append(yyhdyrbts)
                    GSHP_capacity_Y.append(yyhdyrbrl)
                    eff_GSHP_Y.append(yyhdyrbxl)
                    GSHP_minplr_Y.append(yyhdyrbzdfhl)
                    GSHP_pump_ele_Y.append(yyhdyrbsbhdxs)
                    GICE_select_configured_Y.append(yyhrqnrjts)
                    GICE_rated_power_Y.append(yyhrqnrjfdl)
                    effE_GICE_Y.append(yyhrqnrjfdxl)
                    effH_GICE_Y.append(yyhrqnrjzrxl)
                    GICE_minplr_Y.append(yyhrqnrjzdfhl)
                    GB_select_configured_Y.append(yyhrqglts)
                    GB_capacity_Y.append(yyhrqglrl)
                    eff_GB_Y.append(yyhrqglxl)
                    GB_minplr_Y.append(yyhrqglzdfhl)
                    EB_select_configured_Y.append(yyhdglts)
                    EB_capacity_Y.append(yyhdglrl)
                    eff_EB_Y.append(yyhdglxl)
                    EB_minplr_Y.append(yyhdglzdfhl)
                    HC_area_Y.append(yyhgrmj)
                    eff_HC_Y.append(yyhgrxl)
                    WTSH_voulme_Y.append(yyhsxrrl)
                    eff_Hin_Y.append(yyhsxrcnxl)
                    eff_Hout_Y.append(yyhsxrsnxl)
                    WTSH_loss_Y.append(yyhsxrzsl)
                    ratio_Hin_Y.append(yyhsxrcnsl)
                    ratio_Hout_Y.append(yyhsxrsnsl)

                # print(WTSH_voulme_Y)
                # print(eff_Hin_Y)
                # print(eff_Hout_Y)
                # print(WTSH_loss_Y)
                # print(ratio_Hin_Y)
                # print(ratio_Hout_Y)
                get(preselection_regi, preselection_GSHP, preselection_GICE, preselection_GB, preselection_EB,
                    GSHP_select_configured, GSHP_capacity, GSHP_minplr, eff_GSHP, GSHP_pump_ele,
                    GICE_select_configured, GICE_rated_power, effE_GICE, effH_GICE, GICE_minplr,
                    GB_select_configured, GB_capacity, eff_GB, GB_minplr,
                    EB_select_configured, EB_capacity, eff_EB, EB_minplr,
                    preselection_block, preselection_GSHP_X, preselection_GICE_X, preselection_GB_X, preselection_EB_X,
                    preselection_HC_X, preselection_WTSH_X,
                    GSHP_select_configured_X, GSHP_capacity_X, eff_GSHP_X, GSHP_minplr_X, GSHP_pump_ele_X,
                    GICE_select_configured_X, GICE_rated_power_X, effE_GICE_X, effH_GICE_X, GICE_minplr_X,
                    GB_select_configured_X, GB_capacity_X, eff_GB_X, GB_minplr_X,
                    EB_select_configured_X, EB_capacity_X, eff_EB_X, EB_minplr_X,
                    HC_area_X, eff_HC_X,
                    WTSH_voulme_X, eff_Hin_X, eff_Hout_X, WTSH_loss_X, ratio_Hin_X, ratio_Hout_X,
                    eff_pipenet_X, eff_HE_X,
                    preselection_user, preselection_GSHP_Y, preselection_GICE_Y, preselection_GB_Y, preselection_EB_Y,
                    preselection_HC_Y, preselection_WTSH_Y,
                    GSHP_select_configured_Y, GSHP_capacity_Y, eff_GSHP_Y, GSHP_minplr_Y, GSHP_pump_ele_Y,
                    GICE_select_configured_Y, GICE_rated_power_Y, effE_GICE_Y, effH_GICE_Y, GICE_minplr_Y,
                    GB_select_configured_Y, GB_capacity_Y, eff_GB_Y, GB_minplr_Y,
                    EB_select_configured_Y, EB_capacity_Y, eff_EB_Y, EB_minplr_Y,
                    HC_area_Y, eff_HC_Y,
                    WTSH_voulme_Y, eff_Hin_Y, eff_Hout_Y, WTSH_loss_Y, ratio_Hin_Y, ratio_Hout_Y, eff_pipenet_Y,
                    Load_hea, I_solar, Ele_price, Ele_sell_price, Gas_price)

                return True
            else:
                # 排放为优化目标
                return '排放为优化目标'
        else:
            # 未启动
            return '能源站未启动'
    else:
        # 分布式系统
        return '分布式系统'
    conn.close


print(yh('1383', '2020-04-24', '130600', '1'))
