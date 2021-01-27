import matplotlib.pyplot as plt

from granatum_sdk import Granatum
import time

def main():
    tic = time.perf_counter()

    gn = Granatum()

    df = gn.pandas_from_assay(gn.get_import('assay'))
    mingenes = gn.get_arg('min_genes_per_cell')
    maxgenes = gn.get_arg('max_genes_per_cell')
    mt_percent = gn.get_arg('mt_genes_percent')

    gn.add_result("Min unique genes = {}".format(mingenes), "markdown")
    gn.add_result("Max unique genes = {}".format(maxgenes), "markdown")
    gn.add_result("Max mitochondrial percentage = {}".format(mt_percent), "markdown")

    uniquegenecount = df.astype(bool).sum(axis=0)
    totalgenecount = df.sum(axis=0)
    mtrows = df[df.index.str.startswith('MT')]
    mtgenecount = mtrows.sum(axis=0)
    mtpercent = mtgenecount.div(totalgenecount)
    colsmatching = uniquegenecount.T[(uniquegenecount.T >= mingenes) & (uniquegenecount.T <= maxgenes) & (mtpercent.T <= mt_percent)].index.values
    adata = df.loc[:, colsmatching]

    gn.export(gn.assay_from_pandas(adata), "Filtered Cells Assay", dynamic=False)

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished cell filtering step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
