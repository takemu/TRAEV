from config import *
from modeling import *


def prepare_desired_trait_file():
    _, reframed_model = create_gpr_model(model_xml, type='aromatic')
    r_df1 = pd.DataFrame([[r[0], r[1].name] for r in reframed_model.reactions.items()], columns=['r_id','r_name']).set_index('r_id')
    r_df1['DT_PEA'] = 0
    r_df1['DT_BCHA'] = 0
    r_df1.loc[['R_r_1589', 'R_r_2000'], 'DT_PEA'] = 1
    r_df1.loc[['R_r_1580', 'R_r_1581', 'R_r_1862', 'R_r_1865'], 'DT_BCHA'] = 1
    r_df1.to_csv('data/dt_aromatic.csv')

    _, reframed_model = create_gpr_model(model_xml, type='indigoidine')
    r_df2 = pd.DataFrame([[r[0], r[1].name] for r in reframed_model.reactions.items()], columns=['r_id','r_name']).set_index('r_id')
    r_df2['DT_Ind'] = 0
    r_df2.loc['R_r_i003', 'DT_Ind'] = 1
    # r_df2.to_csv('data/dt_indigoidine.csv')

    _, reframed_model = create_gpr_model(model_xml, type='bikaverin')
    r_df3 = pd.DataFrame([[r[0], r[1].name] for r in reframed_model.reactions.items()], columns=['r_id','r_name']).set_index('r_id')
    r_df3['DT_Bik'] = 0
    r_df3.loc['R_r_b003', 'DT_Bik'] = 1
    # r_df3.to_csv('data/dt_bikaverin.csv')

    with pd.ExcelWriter('data/desired_traits.xlsx') as writer:
        r_df1.to_excel(writer, sheet_name='aromatic', index=True)
        r_df2.to_excel(writer, sheet_name='indigoidine', index=True)
        r_df3.to_excel(writer, sheet_name='bikaverin', index=True)


def prepare_ale_environment_file():
    def func(row):
        subs = row['Evolution environment'].split('exchange')
        for i in range(len(subs)):
            subs[i] = subs[i].strip()
            if subs[i] == '':
                del subs[i]
        row['exchange_reactions'] = []
        row['nutrients'] = []
        for i, sub in enumerate(subs):
            name = sub + ' exchange'
            if name == 'glycine exchange':
                name = 'L-glycine exchange'
            row['exchange_reactions'].append(''.join(r_df[r_df['r_name']==name].index.tolist()))
            row['nutrients'].append(name.split(' ')[0])
        row['exchange_reactions'].sort()
        row['nutrients'].sort(key=lambda v: v.upper())
        row['nutrients'] = str(row['nutrients']).strip('[]').replace('\'', '')
        return row
    _, reframed_model = create_gpr_model(model_xml, type='aromatic')
    r_df = pd.DataFrame([[r[0], r[1].name] for r in reframed_model.reactions.items()], columns=['r_id','r_name']).set_index('r_id')
    es_pea = pd.read_excel('data/EvolveX/msb202210980-sup-0001-TableEV1.xlsx', sheet_name='EV_PEA')[['Evolution environment', 'Score']].apply(func, axis=1).rename(columns={'Score': 'EvolveX_score_PEA'})
    es_bcha = pd.read_excel('data/EvolveX/msb202210980-sup-0001-TableEV1.xlsx', sheet_name='EV_BCHA')[['Evolution environment', 'Score']].apply(func, axis=1).rename(columns={'Score': 'EvolveX_score_BCHA'})
    env_df = es_pea.merge(es_bcha[['nutrients', 'EvolveX_score_BCHA']], on='nutrients', how='left').sort_values(by='nutrients')[['nutrients', 'exchange_reactions', 'EvolveX_score_PEA', 'EvolveX_score_BCHA']]
    env_df.index = [f'Env_{i:04d}' for i in range(1, len(env_df) + 1)]
    env_df.index.name = 'env_id'
    env_df.to_csv('data/ale_envs.csv', float_format="%.6f")
    env_df.to_excel('data/ale_envs.xlsx', float_format="%.6f")


if __name__ == '__main__':
    prepare_desired_trait_file()
    prepare_ale_environment_file()
