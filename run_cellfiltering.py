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

    uniquegenecount = df.astype(bool).sum(axis=0)
    totalgenecount = df.sum(axis=0)
    mtrows = df[df.index.str.startswith('MT')]
    mtgenecount = mtrows.sum(axis=0)
    mtpercent = mtgenecount.div(totalgenecount)
    colsmatching = uniquegenecount.T[(uniquegenecount.T >= mingenes) & (uniquegenecount.T <= maxgenes) & (mtpercent.T <= mt_percent)].index.values
    adata = df.loc[:, colsmatching]

    num_gt_min = uniquegenecount.T[(uniquegenecount.T >= mingenes)].index.size
    num_gt_minmax = uniquegenecount.T[(uniquegenecount.T >= mingenes) & (uniquegenecount.T <= maxgenes)].index.size
    num_gt_minmaxmt = uniquegenecount.T[(uniquegenecount.T >= mingenes) & (uniquegenecount.T <= maxgenes) & (mtpercent.T <= mt_percent)].index.size

    gn.add_result("Number of cells is now {} with {} below min genes, and of those, {} above max genes, and of those, {} below mt percentage threshold.".format(len(colsmatching), len(colsmatching) - num_gt_min, num_gt_min - num_gt_minmax, num_gt_minmaxmt - num_gt_minmax), "markdown")

    gn.export(gn.assay_from_pandas(adata), "Filtered Cells Assay", dynamic=False)

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished cell filtering step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
