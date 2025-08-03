# Created by alexandra at 20/12/2023
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.metrics import ConfusionMatrixDisplay, roc_auc_score, roc_curve

def plot_input(data_file):
    data = pd.read_csv(data_file, sep="\t")
    # data["intLabel"] = np.where(data["Label"]=="Cancer", 1, 0)
    data["intStudy"] = np.where(data["study"] == "NeoRheaStudy", 1, data["study"])
    data["intStudy"] = np.where(data["study"] == "PearlStudy", 2, 0)
    X1 = data["short_reads"]
    X2 = data["no_reads"]
    y = data["intStudy"]


    plt.scatter(x=X1,
                y=X2,
                c=y,
                cmap=plt.cm.RdYlBu)

    plt.show()

def draw_result(lst_iter, lst_train, lst_test,title):
    """
    TODO: add comments
    """
    plt.clf()
    plt.plot(lst_iter, lst_train, '-b', label='CrossEntropyLoss Train')
    plt.plot(lst_iter, lst_test, '-r', label='CrossEntropyLoss Test')

    plt.xlabel("Epoch number")
    plt.legend(loc='upper left')
    plt.title(title)
    # save image
    plt.savefig("plots/"+ title+".png")
    plt.clf()


def plot_roc_fpr_tpr(true_y, positive_class_prob, output_file, title, positive_label):
    fpr, tpr, thresholds = roc_curve(true_y, positive_class_prob, pos_label=positive_label)
    lr_auc = np.round(roc_auc_score(true_y, positive_class_prob), 2)
    true_y_np = np.asarray(true_y)
    positive_class_prob_np = np.asarray(positive_class_prob)
    score_confidence_lower, score_confidence_upper = plot_roc_with_ci(true_y_np, positive_class_prob_np, positive_label)

    # plot the roc curve for the model
    gmeans = np.sqrt(tpr * (1 - fpr))
    # locate the index of the largest g-mean
    ix = np.argmax(gmeans)
    plt.plot([0, 1], [0, 1], linestyle='--', label='Chance')
    plt.plot(fpr, tpr, marker='.', label='ALFAssay, AUC=' + str(lr_auc) + "; CI:[" + str(score_confidence_lower) +
                                         "-" + str(score_confidence_upper) + "]")

    plt.scatter(fpr[ix], tpr[ix], marker='o', color='black', label='Best')

    # axis labels
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title, fontsize=7)
    plt.legend(loc=4)
    file = output_file + title + "_roc_fpr_tpr.pdf"
    plt.savefig(file)
    plt.clf()


def plot_roc_fpr_tpr_ichor(true_y, positive_class_prob, output_file, title, positive_label, ichorCNA_pred, caption):

    fpr, tpr, thresholds = roc_curve(true_y, positive_class_prob, pos_label=positive_label)
    lr_auc = np.round(roc_auc_score(true_y, positive_class_prob), 2)
    true_y_np = np.asarray(true_y)
    positive_class_prob_np = np.asarray(positive_class_prob)
    score_confidence_lower, score_confidence_upper = plot_roc_with_ci(true_y_np, positive_class_prob_np, positive_label)

    ichorCNA_pred_bool = np.where(ichorCNA_pred > 0.03, 1, 0)
    fpr_ichorCNA, tpr_ichorCNA, thresholds_ichorCNA = roc_curve(true_y, ichorCNA_pred_bool, pos_label=positive_label)
    lr_auc_ichorCNA = np.round(roc_auc_score(true_y, ichorCNA_pred_bool), 2)

    score_confidence_lower_ichor, score_confidence_upper_ichor = plot_roc_with_ci(true_y, ichorCNA_pred_bool,
                                                                                  positive_label)

    # plot the roc curve for the model
    gmeans = np.sqrt(tpr * (1 - fpr))
    # locate the index of the largest g-mean
    ix = np.argmax(gmeans)
    fig = plt.figure()
    ax1 = fig.add_axes((.1, .2, .8, .7))
    # ax1 = fig.add_axes((1, 1, 1, 1))
    ax1.plot([0, 1], [0, 1], linestyle='--', label='Chance')
    ax1.plot(fpr, tpr, marker='.', label='ALFAssay, AUC=' + str(lr_auc) + "; CI:[" + str(score_confidence_lower) +
             "-" + str(score_confidence_upper) + "]")
    ax1.plot(fpr_ichorCNA, tpr_ichorCNA, marker='.', label='ichorCNA, AUC=' + str(lr_auc_ichorCNA) +
                                                           "; CI:[" + str(score_confidence_lower_ichor) +
                                                           "-" + str(score_confidence_upper_ichor) + "]")

    ax1.scatter(fpr[ix], tpr[ix], marker='o', color='black', label='Best')
    # fig = plt.gcf()
    # axis labels
    # plt.xlabel('''False Positive Rate
    # ''' + caption)
    fig.text(.1, .05, caption)
    plt.ylabel('True Positive Rate')
    plt.title(title, fontsize=7)
    plt.legend(loc=4)
    # fig.text(0.5, -0.05, caption, ha='center')
    file = output_file + title + "_roc_fpr_tpr.pdf"
    plt.savefig(file)
    plt.clf()


def plot_roc_with_ci(true_y, positive_class_prob, pos_label):
    n_bootstraps = 2000
    rng_seed = 42  # control reproducibility
    bootstrapped_scores = []
    bootstrapped_fpr = []
    bootstrapped_tpr = []

    rng = np.random.RandomState(rng_seed)
    for i in range(n_bootstraps):
        # bootstrap by sampling with replacement on the prediction indices
        indices = rng.randint(0, len(positive_class_prob), len(positive_class_prob))
        if len(np.unique(true_y[indices])) < 2:
            # We need at least one positive and one negative sample for ROC AUC
            # to be defined: reject the sample
            continue

        # score = roc_auc_score(y_true[indices], y_pred[indices])
        # fpr, tpr, thresholds = roc_curve(true_y[indices], positive_class_prob[indices], pos_label=pos_label)
        lr_auc = np.round(roc_auc_score(true_y[indices], positive_class_prob[indices]), 2)
        bootstrapped_scores.append(lr_auc)
        # bootstrapped_fpr.append(fpr)
        # bootstrapped_tpr.append(tpr)
        # print("Bootstrap #{} ROC area: {:0.3f}".format(i + 1, lr_auc))
    sorted_scores = np.array(bootstrapped_scores)
    sorted_scores.sort()

    # sorted_bootstrapped_fpr = np.array(bootstrapped_fpr)
    # sorted_bootstrapped_fpr.sort()
    #
    # sorted_bootstrapped_tpr = np.array(bootstrapped_tpr)
    # sorted_bootstrapped_tpr.sort()

    # Computing the lower and upper bound of the 90% confidence interval
    # You can change the bounds percentiles to 0.025 and 0.975 to get
    # a 95% confidence interval instead.
    score_confidence_lower = sorted_scores[int(0.05 * len(sorted_scores))]
    score_confidence_upper = sorted_scores[int(0.95 * len(sorted_scores))]
    # fpr_confidence_lower = sorted_bootstrapped_fpr[int(0.05 * len(sorted_bootstrapped_fpr))]
    # fpr_confidence_upper = sorted_bootstrapped_fpr[int(0.95 * len(sorted_bootstrapped_fpr))]
    # tpr_confidence_lower = sorted_bootstrapped_tpr[int(0.05 * len(sorted_bootstrapped_tpr))]
    # tpr_confidence_upper = sorted_bootstrapped_tpr[int(0.95 * len(sorted_bootstrapped_tpr))]

    # print("Confidence interval for the score: [{:0.3f} - {:0.3}]".format(
    #     score_confidence_lower, score_confidence_upper))

    return score_confidence_lower, score_confidence_upper #, \
           # fpr_confidence_lower, fpr_confidence_upper, \
           # tpr_confidence_lower, tpr_confidence_upper





def plot_roc_fpr_tpr_multimodels(true_y, ichor_positive_class_prob,
                                 fragle_positive_class_prob,
                                 alfassay_positive_class_prob,
                                 alfassay_ichor_positive_class_prob,
                                 output_file, title, positive_label, caption):
    #alfassay
    fpr_alfassay, tpr_alfassay, thresholds_alfassay = roc_curve(true_y, alfassay_positive_class_prob, pos_label=positive_label)
    lr_auc_alfassay= np.round(roc_auc_score(true_y, alfassay_positive_class_prob), 2)
    true_y_np = np.asarray(true_y)
    positive_class_prob_np = np.asarray(alfassay_positive_class_prob)
    score_confidence_lower_alfassay, score_confidence_upper_alfassay = plot_roc_with_ci(true_y_np, positive_class_prob_np, positive_label)

    #ichor
    fpr_ichorCNA, tpr_ichorCNA, thresholds_ichorCNA = roc_curve(true_y, ichor_positive_class_prob, pos_label=positive_label)
    lr_auc_ichorCNA = np.round(roc_auc_score(true_y, ichor_positive_class_prob), 2)

    score_confidence_lower_ichor, score_confidence_upper_ichor = plot_roc_with_ci(true_y, ichor_positive_class_prob,
                                                                                  positive_label)

    #fragle
    fpr_fragle, tpr_fragle, thresholds_fragle = roc_curve(true_y, fragle_positive_class_prob,
                                                                pos_label=positive_label)
    lr_auc_fragle = np.round(roc_auc_score(true_y, fragle_positive_class_prob), 2)

    score_confidence_lower_fragle, score_confidence_upper_fragle = plot_roc_with_ci(true_y, fragle_positive_class_prob,
                                                                                  positive_label)

    # ichor + alfassay
    fpr_ichor_alfa, tpr_ichor_alfa, thresholds_ichor_alfa = roc_curve(true_y, alfassay_ichor_positive_class_prob,
                                                          pos_label=positive_label)
    lr_auc_ichor_alfa = np.round(roc_auc_score(true_y, alfassay_ichor_positive_class_prob), 2)

    score_confidence_lower_ichor_alfa, score_confidence_upper_ichor_alfa = plot_roc_with_ci(true_y, alfassay_ichor_positive_class_prob,
                                                                                    positive_label)

    # plot the roc curve for the model
    gmeans = np.sqrt(tpr_alfassay * (1 - fpr_alfassay))
    # locate the index of the largest g-mean
    ix = np.argmax(gmeans)
    fig = plt.figure()
    ax1 = fig.add_axes((.1, .2, .8, .7))
    # ax1 = fig.add_axes((1, 1, 1, 1))
    ax1.plot([0, 1], [0, 1], linestyle='--', label='Chance')
    ax1.plot(fpr_alfassay, tpr_alfassay, marker='.', label='ALFAssay, AUC=' + str(lr_auc_alfassay) + "; CI:[" + str(score_confidence_lower_alfassay) +
             "-" + str(score_confidence_upper_alfassay) + "]")
    ax1.plot(fpr_ichorCNA, tpr_ichorCNA, marker='.', label='ichorCNA, AUC=' + str(lr_auc_ichorCNA) +
                                                           "; CI:[" + str(score_confidence_lower_ichor) +
                                                           "-" + str(score_confidence_upper_ichor) + "]")

    ax1.plot(fpr_fragle, tpr_fragle, marker='.', label='Fragle, AUC=' + str(lr_auc_fragle) +
                                                           "; CI:[" + str(score_confidence_lower_fragle) +
                                                           "-" + str(score_confidence_upper_fragle) + "]")

    ax1.plot(fpr_ichor_alfa, tpr_ichor_alfa, marker='.', label='ALFAssay + ichorCNA, AUC=' + str(lr_auc_ichor_alfa) +
                                                       "; CI:[" + str(score_confidence_lower_ichor_alfa) +
                                                       "-" + str(score_confidence_upper_ichor_alfa) + "]")

    ax1.scatter(fpr_alfassay[ix], tpr_alfassay[ix], marker='o', color='black', label='Best')
    # fig = plt.gcf()
    # axis labels
    # plt.xlabel('''False Positive Rate
    # ''' + caption)
    fig.text(.1, .05, caption)
    plt.ylabel('True Positive Rate')
    plt.title(title, fontsize=7)
    plt.legend(loc=4)
    # fig.text(0.5, -0.05, caption, ha='center')
    file = output_file + title + "_roc_fpr_tpr.pdf"
    plt.savefig(file)
    plt.clf()